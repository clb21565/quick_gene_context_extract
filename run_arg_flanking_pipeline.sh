#!/usr/bin/env bash
# run_arg_flanking_pipeline.sh
#
# ARG identification + flanking-region extraction pipeline:
#   contigs -> Prodigal ORFs -> DIAMOND vs ARG database -> pad N kb each
#   side of every hit -> extract sequence -> one FASTA per ARG.
#
# See README.md for full documentation.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage: run_arg_flanking_pipeline.sh -a ASSEMBLY.fasta -d ARG_DB.dmnd -o OUTDIR [options]

Required:
  -a, --assembly FILE        Nucleotide assembly/contigs FASTA
  -d, --arg-db FILE          DIAMOND-formatted ARG protein database (.dmnd)
                              (not required if --skip-diamond is set)

Common options:
  -o, --outdir DIR            Output directory (default: arg_flank_results)
  -p, --prefix NAME           Sample/output prefix (default: sample)
  -f, --flank-kb N            Flank size in kb to extract on each side of every ARG hit (default: 5)
  -m, --merge-distance N      Merge ARG hits within N bp of each other before padding (default: 0, no merging)
  -e, --evalue N               DIAMOND e-value threshold (default: 1e-10)
  -T, --threads N              Threads for DIAMOND (default: 4)

ARG-label naming for split FASTA files (optional, see README):
  --stitle-sep CHAR             Field separator in DIAMOND stitle (e.g. '|' for CARD-style headers)
  --stitle-field N               1-based field index to use as the ARG label after splitting by --stitle-sep

Skip steps if you already have outputs from elsewhere in your workflow:
  --skip-prodigal                Skip Prodigal; requires --orfs
  --orfs FILE                     Existing Prodigal protein FASTA (prodigal -p meta headers)
  --skip-diamond                  Skip DIAMOND; requires --diamond-tsv
  --diamond-tsv FILE                Existing DIAMOND outfmt6 TSV (qtitle stitle pident bitscore evalue, no header)

  -h, --help                      Show this help
EOF
}

# ---- defaults ----
OUTDIR="arg_flank_results"
PREFIX="sample"
FLANK_KB=5
MERGE_DIST=0
EVALUE="1e-10"
THREADS=4
SKIP_PRODIGAL=0
SKIP_DIAMOND=0
ORF_FAA=""
DIAMOND_TSV=""
STITLE_SEP=""
STITLE_FIELD=""
ASSEMBLY=""
ARG_DB=""

# ---- argument parsing ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    -a|--assembly) ASSEMBLY="$2"; shift 2 ;;
    -d|--arg-db) ARG_DB="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    -p|--prefix) PREFIX="$2"; shift 2 ;;
    -f|--flank-kb) FLANK_KB="$2"; shift 2 ;;
    -m|--merge-distance) MERGE_DIST="$2"; shift 2 ;;
    -e|--evalue) EVALUE="$2"; shift 2 ;;
    -T|--threads) THREADS="$2"; shift 2 ;;
    --stitle-sep) STITLE_SEP="$2"; shift 2 ;;
    --stitle-field) STITLE_FIELD="$2"; shift 2 ;;
    --skip-prodigal) SKIP_PRODIGAL=1; shift ;;
    --orfs) ORF_FAA="$2"; shift 2 ;;
    --skip-diamond) SKIP_DIAMOND=1; shift ;;
    --diamond-tsv) DIAMOND_TSV="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

# ---- validation ----
[[ -z "$ASSEMBLY" ]] && { echo "ERROR: --assembly is required" >&2; exit 1; }
[[ ! -f "$ASSEMBLY" ]] && { echo "ERROR: assembly file not found: $ASSEMBLY" >&2; exit 1; }
if [[ "$SKIP_DIAMOND" -eq 0 && -z "$ARG_DB" ]]; then
  echo "ERROR: --arg-db is required unless --skip-diamond is set" >&2; exit 1
fi
if [[ "$SKIP_DIAMOND" -eq 1 && -z "$DIAMOND_TSV" ]]; then
  echo "ERROR: --diamond-tsv is required when --skip-diamond is set" >&2; exit 1
fi
if [[ "$SKIP_PRODIGAL" -eq 1 && -z "$ORF_FAA" ]]; then
  echo "ERROR: --orfs is required when --skip-prodigal is set" >&2; exit 1
fi

ASSEMBLY="$(realpath "$ASSEMBLY")"
[[ -n "$ARG_DB" ]] && ARG_DB="$(realpath "$ARG_DB")"
[[ -n "$ORF_FAA" ]] && ORF_FAA="$(realpath "$ORF_FAA")"
[[ -n "$DIAMOND_TSV" ]] && DIAMOND_TSV="$(realpath "$DIAMOND_TSV")"

mkdir -p "$OUTDIR"
OUTDIR="$(realpath "$OUTDIR")"
LOG="$OUTDIR/${PREFIX}.pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "=== ARG flanking-region extraction pipeline ==="
echo "Assembly:        $ASSEMBLY"
echo "Output dir:      $OUTDIR"
echo "Flank size:      ${FLANK_KB} kb each side"
echo "Merge distance:  ${MERGE_DIST} bp"
echo

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: required tool not found in PATH: $1" >&2; exit 1; }; }
need bedtools
need samtools
need python3
[[ "$SKIP_PRODIGAL" -eq 0 ]] && need prodigal
[[ "$SKIP_DIAMOND" -eq 0 ]] && need diamond

# faidx the assembly via a symlink in OUTDIR, in case the source dir is read-only
ASSEMBLY_LINK="$OUTDIR/$(basename "$ASSEMBLY")"
[[ -e "$ASSEMBLY_LINK" ]] || ln -s "$ASSEMBLY" "$ASSEMBLY_LINK"
samtools faidx "$ASSEMBLY_LINK"
cut -f1,2 "${ASSEMBLY_LINK}.fai" > "$OUTDIR/${PREFIX}.genome"

# --- Step 1: Prodigal ORF calling ---
if [[ "$SKIP_PRODIGAL" -eq 0 ]]; then
  echo "[1/7] Calling ORFs with Prodigal..."
  ORF_FAA="$OUTDIR/${PREFIX}.faa"
  prodigal -i "$ASSEMBLY" -a "$ORF_FAA" -p meta -q
else
  echo "[1/7] Skipping Prodigal, using existing ORFs: $ORF_FAA"
fi

# --- Step 2: DIAMOND search against ARG database ---
if [[ "$SKIP_DIAMOND" -eq 0 ]]; then
  echo "[2/7] Searching ORFs against ARG database with DIAMOND..."
  DIAMOND_TSV="$OUTDIR/${PREFIX}.diamond.tsv"
  diamond blastp \
    -q "$ORF_FAA" \
    -d "$ARG_DB" \
    -o "$DIAMOND_TSV" \
    --outfmt 6 qtitle stitle pident bitscore evalue \
    --evalue "$EVALUE" \
    --threads "$THREADS" \
    --max-target-seqs 1
else
  echo "[2/7] Skipping DIAMOND, using existing hits: $DIAMOND_TSV"
fi

if [[ ! -s "$DIAMOND_TSV" ]]; then
  echo "No ARG hits found (empty DIAMOND output: $DIAMOND_TSV). Stopping."
  exit 0
fi

# --- Step 3: parse hits -> processed.tsv + BED ---
echo "[3/7] Parsing DIAMOND hits and Prodigal coordinates..."
python3 "$SCRIPT_DIR/parse_diamond_prodigal.py" \
  -i "$DIAMOND_TSV" -p "$PREFIX" -o "$OUTDIR"

PROCESSED_TSV="$OUTDIR/${PREFIX}.processed.tsv"
BED="$OUTDIR/${PREFIX}.bed"

# --- Step 4: sort, optionally merge nearby hits ---
echo "[4/7] Sorting BED file..."
sort -k1,1 -k2,2n "$BED" > "$OUTDIR/${PREFIX}.sorted.bed"

if [[ "$MERGE_DIST" -gt 0 ]]; then
  echo "      Merging ARG hits within ${MERGE_DIST} bp of each other..."
  bedtools merge -i "$OUTDIR/${PREFIX}.sorted.bed" -d "$MERGE_DIST" -c 4 -o collapse \
    > "$OUTDIR/${PREFIX}.merged_regions.bed"
  REGIONS_BED="$OUTDIR/${PREFIX}.merged_regions.bed"
else
  REGIONS_BED="$OUTDIR/${PREFIX}.sorted.bed"
fi

# --- Step 5: pad regions by flank size ---
echo "[5/7] Padding regions by ${FLANK_KB} kb on each side..."
FLANK_BP=$(( FLANK_KB * 1000 ))
bedtools slop -i "$REGIONS_BED" -g "$OUTDIR/${PREFIX}.genome" -b "$FLANK_BP" \
  > "$OUTDIR/${PREFIX}.padded_regions.bed"

# --- Step 6: extract padded region sequences ---
echo "[6/7] Extracting flanking sequences..."
bedtools getfasta -fi "$ASSEMBLY_LINK" -bed "$OUTDIR/${PREFIX}.padded_regions.bed" -name \
  -fo "$OUTDIR/${PREFIX}.extracted_regions.fasta"

# --- Step 7: split into one FASTA per ARG hit ---
echo "[7/7] Splitting into per-ARG FASTA files..."
SPLIT_ARGS=()
[[ -n "$STITLE_SEP" ]] && SPLIT_ARGS+=(--stitle-sep "$STITLE_SEP")
[[ -n "$STITLE_FIELD" ]] && SPLIT_ARGS+=(--stitle-field "$STITLE_FIELD")

python3 "$SCRIPT_DIR/split_by_arg.py" \
  -f "$OUTDIR/${PREFIX}.extracted_regions.fasta" \
  -t "$PROCESSED_TSV" \
  -o "$OUTDIR/split_arg_fastas" \
  "${SPLIT_ARGS[@]}"

echo
echo "✅ Done. Key outputs in $OUTDIR/:"
echo "   ${PREFIX}.processed.tsv           - ARG hit table with ORF coordinates"
echo "   ${PREFIX}.padded_regions.bed      - flanking regions (padded +/- ${FLANK_KB}kb)"
echo "   ${PREFIX}.extracted_regions.fasta - flanking sequences, one record per region"
echo "   split_arg_fastas/                 - one FASTA per ARG, named from the DIAMOND hit"
