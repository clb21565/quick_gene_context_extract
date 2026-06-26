#!/usr/bin/env python3
"""
parse_diamond_prodigal.py

Parse DIAMOND tabular hits (outfmt 6, with full qtitle/stitle headers) where
the queries are proteins called by Prodigal in meta mode (`prodigal -p meta`),
and produce:

  1) <prefix>.processed.tsv - one row per ARG hit, with ORF coordinates pulled
     out of the Prodigal header (contig, start, stop, strand, partial,
     start_codon, rbs_motif, rbs_spacer, gc_cont) plus the DIAMOND hit info.
  2) <prefix>.bed - BED6 of each ORF hit's location on the original contig.

Prodigal protein FASTA headers look like:
  contig_1_1 # 61 # 423 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.502

Fields are separated by " # ", which is what this script splits on (this is
more robust than regex-matching across the header).

Usage:
  parse_diamond_prodigal.py -i diamond_hits.tsv -p sample_prefix -o outdir
"""
import argparse
import os
import sys

import pandas as pd

COLUMNS = ["qtitle", "stitle", "pident", "bitscore", "evalue"]


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "-i", "--diamond-tsv", required=True,
        help="DIAMOND outfmt 6 TSV with columns: qtitle stitle pident bitscore evalue (no header row)",
    )
    ap.add_argument("-p", "--prefix", default="sample", help="Output file prefix")
    ap.add_argument("-o", "--outdir", default=".", help="Output directory")
    return ap.parse_args()


def parse_prodigal_header(qtitle):
    """Split a Prodigal protein-FASTA header into coordinates + metadata."""
    parts = qtitle.split(" # ")
    if len(parts) < 4:
        raise ValueError(f"qtitle does not look like a Prodigal header: {qtitle!r}")
    orf_id = parts[0].split(" ")[0]
    start = int(parts[1])
    stop = int(parts[2])
    strand = "+" if parts[3].strip() == "1" else "-"
    meta = {}
    if len(parts) > 4:
        for kv in parts[4].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                meta[k.strip()] = v.strip()
    # Prodigal names ORFs "<contig>_<orf_number>" -- contig is everything
    # before the final underscore-delimited field.
    contig = orf_id.rsplit("_", 1)[0]
    return contig, orf_id, start, stop, strand, meta


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.diamond_tsv, sep="\t", header=None, names=COLUMNS)
    if df.empty:
        sys.exit(f"No rows found in {args.diamond_tsv} -- no ARG hits to process.")

    parsed = df["qtitle"].apply(parse_prodigal_header)
    df["contig"] = [p[0] for p in parsed]
    df["orf"] = [p[1] for p in parsed]
    df["start"] = [p[2] for p in parsed]
    df["stop"] = [p[3] for p in parsed]
    df["strand"] = [p[4] for p in parsed]
    df["partial"] = [p[5].get("partial", "") for p in parsed]
    df["start_codon"] = [p[5].get("start_type", "") for p in parsed]
    df["rbs_motif"] = [p[5].get("rbs_motif", "") for p in parsed]
    df["rbs_spacer"] = [p[5].get("rbs_spacer", "") for p in parsed]
    df["gc_cont"] = [p[5].get("gc_cont", "") for p in parsed]

    out_cols = [
        "contig", "orf", "start", "stop", "strand",
        "stitle", "pident", "bitscore", "evalue",
        "partial", "start_codon", "rbs_motif", "rbs_spacer", "gc_cont",
    ]
    processed_path = os.path.join(args.outdir, f"{args.prefix}.processed.tsv")
    df[out_cols].to_csv(processed_path, sep="\t", index=False)

    bed_path = os.path.join(args.outdir, f"{args.prefix}.bed")
    with open(bed_path, "w") as bed:
        for _, row in df.iterrows():
            # BED is 0-based, half-open
            bed.write(
                f"{row['contig']}\t{row['start'] - 1}\t{row['stop']}\t"
                f"{row['orf']}\t0\t{row['strand']}\n"
            )

    print(f"Wrote {len(df)} ARG hits -> {processed_path}, {bed_path}")


if __name__ == "__main__":
    main()
