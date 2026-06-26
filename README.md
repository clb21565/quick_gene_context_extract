# ARG Flanking-Region Extraction Pipeline

A small pipeline that, starting from an assembly, finds antibiotic resistance
genes (ARGs) and extracts the genomic neighborhood around each one (N kb to
the left and right), writing out one FASTA per ARG. Useful for downstream
work on mobile genetic element context, co-localization with MGEs/MGE
markers, synteny comparisons, etc.

## Pipeline overview

```
assembly.fasta
   |
   v
[1] Prodigal (-p meta)         -> ORF protein FASTA (with coordinate headers)
   |
   v
[2] DIAMOND blastp vs ARG DB   -> ARG hit table (qtitle, stitle, pident, bitscore, evalue, qcovhsp)
   |     (filtered by --min-pident / --min-qcov / --evalue)
   |
   v
[3] parse_diamond_prodigal.py  -> processed.tsv (hit table + ORF coords) + BED
   |
   v
[4] sort / (optional) bedtools merge nearby hits
   |
   v
[5] bedtools slop               -> pad each region +/- N kb (clipped at contig ends)
   |
   v
[6] bedtools getfasta           -> one FASTA record per padded region
   |
   v
[7] split_by_arg.py             -> split_arg_fastas/<ARG_label>.fasta, one file per ARG
```

## Files

| File | Purpose |
|---|---|
| `run_arg_flanking_pipeline.sh` | Main driver. Runs all steps end to end. |
| `parse_diamond_prodigal.py` | Parses DIAMOND hits + Prodigal headers into a tidy TSV and BED. |
| `split_by_arg.py` | Splits the extracted multi-record FASTA into one file per ARG. |

## Requirements

- [Prodigal](https://github.com/hyattpd/Prodigal) (skip with `--skip-prodigal` if you already call ORFs elsewhere)
- [DIAMOND](https://github.com/bbuchfink/diamond) (skip with `--skip-diamond` if you already have hits)
- [BEDTools](https://bedtools.readthedocs.io/) >= 2.27 (for the `name::chrom:start-stop` FASTA header format used by `getfasta -name`; older versions also work, see note below)
- [samtools](http://www.htslib.org/)
- Python 3 with `pandas`

On TinkerCliffs or similar HPC, build the env from `environment.yml`:

```bash
conda env create -f environment.yml
conda activate arg-flank
```

## Quick start

```bash
./run_arg_flanking_pipeline.sh \
  -a assembly.fasta \
  -d card_protein_homolog.dmnd \
  -o results/ \
  -p mysample \
  -f 5 \
  -m 0
```

This calls ORFs with Prodigal, searches them against the CARD database with
DIAMOND, then extracts 5 kb of flanking sequence on each side of every ARG
hit into `results/split_arg_fastas/<label>.fasta`.

### If you already have ORFs and/or DIAMOND hits

```bash
./run_arg_flanking_pipeline.sh \
  -a assembly.fasta \
  --skip-prodigal --orfs proteins.faa \
  --skip-diamond --diamond-tsv hits.tsv \
  -o results/ -p mysample -f 10
```

`--orfs` must still be a Prodigal-style protein FASTA (used only as the
source of headers when DIAMOND was run separately), and `--diamond-tsv`
must have columns `qtitle stitle pident bitscore evalue`, optionally
followed by `qcovhsp` (6 columns total) -- i.e. DIAMOND was run with:

```bash
diamond blastp -q proteins.faa -d arg_db.dmnd -o hits.tsv \
  --outfmt 6 qtitle stitle pident bitscore evalue qcovhsp \
  --id 80 --query-cover 80
```

The parser auto-detects whether `qcovhsp` is present (5 vs. 6 columns), so
older 5-column TSVs from `--skip-diamond` still work -- they just won't have
a `qcovhsp` column in `processed.tsv`.

`qtitle`/`stitle` (rather than the default `qseqid`/`sseqid`) are required
because the full Prodigal header -- including the `# start # stop # strand #
ID=...` coordinate block -- has to survive into the DIAMOND output for the
coordinate parsing step to work.

## Arguments

```
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
  -I, --min-pident N            Minimum percent identity for a DIAMOND hit (default: 80)
  -Q, --min-qcov N               Minimum percent query coverage for a DIAMOND hit (default: 80)
  -T, --threads N              Threads for DIAMOND (default: 4)

ARG-label naming for split FASTA files:
  --stitle-sep CHAR             Field separator in DIAMOND stitle (e.g. '|' for CARD-style headers)
  --stitle-field N               1-based field index to use as the ARG label after splitting by --stitle-sep

Skip steps if you already have outputs from elsewhere in your workflow:
  --skip-prodigal                Skip Prodigal; requires --orfs
  --orfs FILE                     Existing Prodigal protein FASTA (prodigal -p meta headers)
  --skip-diamond                  Skip DIAMOND; requires --diamond-tsv
  --diamond-tsv FILE                Existing DIAMOND outfmt6 TSV (qtitle stitle pident bitscore evalue, no header)
```

## Output

```
<outdir>/
├── <prefix>.pipeline.log              # full run log
├── <prefix>.faa                       # Prodigal protein FASTA (unless --skip-prodigal)
├── <prefix>.diamond.tsv               # raw DIAMOND hits (unless --skip-diamond)
├── <prefix>.processed.tsv             # one row per ARG hit: contig, ORF coords, strand, hit stats
├── <prefix>.bed                       # per-hit BED (unsorted)
├── <prefix>.sorted.bed                # sorted BED
├── <prefix>.merged_regions.bed        # only if -m/--merge-distance > 0
├── <prefix>.genome                    # contig lengths, for bedtools slop/getfasta
├── <prefix>.padded_regions.bed        # regions after +/- N kb padding (clipped at contig ends)
├── <prefix>.extracted_regions.fasta   # one record per padded region
└── split_arg_fastas/
    └── <ARG_label>.fasta              # one FASTA per ARG, containing its flanking region(s)
```

## Notes / customization

- **Merging nearby hits (`-m`).** With the default `-m 0`, each ARG hit gets
  its own independent +/- N kb window. If you set `-m` to a nonzero distance,
  ARG hits within that distance of each other on the same contig are merged
  into a single region before padding (useful for clusters of resistance
  genes, e.g. on an integron or AMR cassette). `split_by_arg.py` handles
  both cases the same way: if a merged region contains more than one ARG,
  its sequence is written into each of those ARGs' output FASTA files.

- **Contig edges.** `bedtools slop -g <prefix>.genome` automatically clips
  padding at contig boundaries instead of producing invalid coordinates, so
  short contigs or ARGs near a contig end just get whatever flanking
  sequence is actually available -- no wraparound (treated as linear, not
  circular).

- **Naming split FASTA files (`--stitle-sep` / `--stitle-field`).** By
  default the ARG label is just the first whitespace-delimited token of the
  DIAMOND `stitle`. Many ARG databases use pipe-delimited headers, e.g. CARD:

  ```
  gb|AAA12345|ARO:3000001|mecA
  ```

  Use `--stitle-sep '|' --stitle-field 2` to name files by accession
  (`AAA12345.fasta`) or `--stitle-field 4` to name them by gene symbol
  (`mecA.fasta`). Adjust the separator/field index to match whatever
  database you're searching against.

- **Identity/coverage cutoffs (`-I`/`-Q`).** Defaults are 80% identity and
  80% query coverage -- a common starting point for ARG screening tools
  (e.g. ABRicate), meant to avoid flanking-region extraction around marginal
  or partial-length hits. These are passed straight to DIAMOND's `--id` and
  `--query-cover`, so filtering happens at search time, not after the fact.
  Tighten or loosen them for your database; CARD's curated protein homolog
  models, for instance, often warrant a stricter cutoff than the loosely
  related "variant" models in the same database.

- **Multiple DIAMOND hits per ORF.** The driver script runs DIAMOND with
  `--max-target-seqs 1` (best hit only). If you want to keep multiple hits
  per ORF (e.g. to flag ambiguous/multi-domain calls), edit that flag in
  `run_arg_flanking_pipeline.sh`; both downstream Python scripts already
  handle multiple rows per ORF correctly.

- **Older BEDTools.** Versions before ~2.27 don't append `::chrom:start-stop`
  to `-name` headers. `split_by_arg.py` handles this fine either way since it
  only looks at everything before the first `::`, falling back to the whole
  header if there isn't one.
