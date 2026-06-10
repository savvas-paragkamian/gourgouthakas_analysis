#!/usr/bin/env python3
"""isoQC replica - quality-trim Sanger .ab1 traces (isolateR::isoQC).

Two trimming modes:
  auto  (default)  Mott's algorithm (Biopython "abi-trim"), the isolateR
                   "auto cutoff (recommended)" behaviour.
  window           fixed sliding window: trim ends while the mean Phred over
                   --window bases stays below --phred (isolateR defaults Q20/15).

Reads from all *.ab1 in the input dir, writes a trimmed multi-FASTA and a
per-read QC table. Reads shorter than --min-len after trimming are FAIL.

Usage:
  01_isoqc.py --in <ab1_dir> --fasta <out.fasta> --csv <out.csv> \
      [--mode auto|window] [--phred 20] [--window 15] [--min-len 200]
"""
import argparse
import csv
import glob
import os
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def mean_q(quals):
    return sum(quals) / len(quals) if quals else 0.0


def trim_window(seq, quals, phred, window):
    """Trim both ends while the windowed mean Phred is below `phred`."""
    n = len(quals)
    if n < window:
        return 0, n
    # left
    start = 0
    while start <= n - window and mean_q(quals[start:start + window]) < phred:
        start += 1
    # right
    end = n
    while end >= window and mean_q(quals[end - window:end]) < phred:
        end -= 1
    if end <= start:
        return 0, 0
    return start, end


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--mode", choices=["auto", "window"], default="auto")
    ap.add_argument("--phred", type=float, default=20.0)
    ap.add_argument("--window", type=int, default=15)
    ap.add_argument("--min-len", type=int, default=200)
    args = ap.parse_args()

    ab1s = sorted(glob.glob(os.path.join(args.indir, "*.ab1")))
    if not ab1s:
        sys.exit(f"No .ab1 files in {args.indir}")

    rows = []
    n_pass = 0
    with open(args.fasta, "w") as fa:
        for path in ab1s:
            rid = os.path.splitext(os.path.basename(path))[0]
            try:
                rec = SeqIO.read(path, "abi")
            except Exception as e:  # noqa: BLE001 - report and skip bad traces
                rows.append([rid, 0, 0, 0.0, "FAIL", f"read_error:{e}"])
                continue

            quals = rec.letter_annotations.get("phred_quality", [])
            raw_len = len(rec.seq)

            if args.mode == "auto":
                trec = SeqIO.read(path, "abi-trim")  # Mott's algorithm
                tseq = str(trec.seq)
                # recover the trimmed-region qualities for reporting
                tq = trec.letter_annotations.get("phred_quality", [])
            else:
                s, e = trim_window(str(rec.seq), quals, args.phred, args.window)
                tseq = str(rec.seq)[s:e]
                tq = quals[s:e]

            tlen = len(tseq)
            mq = round(mean_q(tq), 2)
            if tlen >= args.min_len:
                decision, note = "PASS", ""
                SeqIO.write(
                    SeqRecord(trec.seq if args.mode == "auto" else rec.seq[s:e],
                              id=rid, description=""),
                    fa, "fasta")
                n_pass += 1
            else:
                decision, note = "FAIL", f"len<{args.min_len}"
            rows.append([rid, raw_len, tlen, mq, decision, note])

    with open(args.csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "raw_len", "trim_len", "mean_phred",
                    "decision", "note"])
        w.writerows(rows)

    print(f"[isoQC:{args.mode}] {len(ab1s)} reads -> {n_pass} PASS "
          f"-> {args.fasta}")


if __name__ == "__main__":
    main()
