#!/usr/bin/env python3
"""
split_by_arg.py

Split a multi-record FASTA of extracted flanking regions (produced by
`bedtools getfasta -name`) into one FASTA file per ARG, named using the
DIAMOND hit's `stitle` field.

Each input record's header is the BED "name" field carried through the
pipeline: a comma-separated list of Prodigal ORF IDs covered by that region
(more than one if regions were merged with `bedtools merge -c 4 -o collapse`),
optionally followed by "::contig:start-stop" (bedtools >= 2.27 appends this
automatically when using `-name`).

If a region covers more than one distinct ARG (only possible when hits were
merged), its sequence is written to every relevant per-ARG FASTA file, with
a note of which other ARGs share that window.

Usage:
  split_by_arg.py -f extracted_regions.fasta -t sample.processed.tsv \
      -o split_arg_fastas/ [--stitle-sep '|'] [--stitle-field 2]

For CARD-style headers ("gb|ACCESSION|ARO:3002999|geneName ..."), use
  --stitle-sep '|' --stitle-field 2   (accession)
  --stitle-sep '|' --stitle-field 4   (gene name)
"""
import argparse
import os
import re
import sys

import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-f", "--fasta", required=True,
                     help="Extracted regions FASTA (bedtools getfasta -name output)")
    ap.add_argument("-t", "--processed-tsv", required=True,
                     help="<prefix>.processed.tsv produced by parse_diamond_prodigal.py")
    ap.add_argument("-o", "--outdir", default="split_arg_fastas",
                     help="Output directory for per-ARG FASTA files")
    ap.add_argument("--stitle-sep", default=None,
                     help="Field separator in the DIAMOND stitle (e.g. '|' for CARD-style headers). "
                          "If omitted, the first whitespace-delimited token of stitle is used as the label.")
    ap.add_argument("--stitle-field", type=int, default=None,
                     help="1-based field index to use after splitting stitle by --stitle-sep")
    return ap.parse_args()


def sanitize(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name.strip())


def read_fasta(path):
    """Minimal FASTA reader; yields (header, sequence)."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def build_orf_to_label(tsv_path, stitle_sep, stitle_field):
    df = pd.read_csv(tsv_path, sep="\t")
    labels = {}
    for _, row in df.iterrows():
        stitle = str(row["stitle"])
        if stitle_sep and stitle_field:
            fields = stitle.split(stitle_sep)
            idx = stitle_field - 1
            label = fields[idx] if 0 <= idx < len(fields) else stitle
        else:
            toks = stitle.split()
            label = toks[0] if toks else stitle
        labels[row["orf"]] = sanitize(label)
    return labels


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    orf_to_label = build_orf_to_label(args.processed_tsv, args.stitle_sep, args.stitle_field)

    handles = {}
    n_records = 0
    n_written = 0
    try:
        for header, seq in read_fasta(args.fasta):
            n_records += 1
            name_field = header.split("::")[0]
            orf_ids = name_field.split(",")
            labels = sorted({orf_to_label[o] for o in orf_ids if o in orf_to_label})
            if not labels:
                print(f"WARNING: no ARG label found for region '{header}', skipping", file=sys.stderr)
                continue
            for label in labels:
                if label not in handles:
                    handles[label] = open(os.path.join(args.outdir, f"{label}.fasta"), "w")
                handles[label].write(f">{header}\n{seq}\n")
                n_written += 1
    finally:
        for h in handles.values():
            h.close()

    print(f"Processed {n_records} extracted regions -> wrote {n_written} sequences "
          f"across {len(handles)} ARG FASTA files in {args.outdir}/")


if __name__ == "__main__":
    main()
