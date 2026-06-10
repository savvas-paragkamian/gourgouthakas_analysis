#!/usr/bin/env python3
"""Merge the SILVA, NCBI and GTDB taxonomy of every isolate into one TSV.

Reads the per-batch isoTAX outputs (results/sanger/<batch>/<batch>.tax.<db>.csv)
rather than the merged ones: read ids restart per plate (e.g. 1_27f-A01 exists
on several plates) and the merged step keeps only the first of each collision,
so the per-batch files are the complete, unambiguous source. Each microbe is
keyed by (plate, microbe_id); for every database the reference identifier
(best hit), % identity and lineage are reported side by side.

Usage:
  08_taxonomy_table.py [--results-dir results/sanger] \
                       [--out results/taxonomy_per_microbe.tsv]
"""
import argparse
import csv
import glob
import os

DBS = ["silva", "ncbi", "gtdb"]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def read_tax(path):
    """query -> dict(best_hit, pct_id, genus, species, lineage)."""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            lineage = ";".join(row[r] for r in RANKS if row.get(r))
            out[row["query"]] = {
                "id": row.get("best_hit", ""),
                "pct": row.get("pct_id", ""),
                "genus": row.get("genus", ""),
                "species": row.get("species", ""),
                "lineage": lineage,
            }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results/sanger")
    ap.add_argument("--out", default="results/taxonomy_per_microbe.tsv")
    ap.add_argument("--acc2taxid", default="data/ref/ncbi_16S.acc2taxid.tsv",
                    help="NCBI accession<TAB>taxid map (from setup_db.sh)")
    args = ap.parse_args()

    # NCBI accession -> taxid (the best_hit in *.tax.ncbi.csv is the accession)
    acc2taxid = {}
    if os.path.isfile(args.acc2taxid):
        with open(args.acc2taxid) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    acc2taxid[p[0]] = p[1]

    batches = sorted(
        d for d in glob.glob(os.path.join(args.results_dir, "*"))
        if os.path.isdir(d) and os.path.basename(d) != "merged")

    header = ["plate", "microbe_id"]
    for db in DBS:
        header.append(f"{db}_id")
        if db == "ncbi":
            header.append("ncbi_taxid")
        header += [f"{db}_pct_id", f"{db}_genus", f"{db}_species",
                   f"{db}_lineage"]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    n_rows = 0
    with open(args.out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(header)
        for bdir in batches:
            plate = os.path.basename(bdir)
            tax = {}
            for db in DBS:
                p = os.path.join(bdir, f"{plate}.tax.{db}.csv")
                tax[db] = read_tax(p) if os.path.isfile(p) else {}
            microbes = sorted(set().union(*(t.keys() for t in tax.values())))
            for mid in microbes:
                row = [plate, mid]
                for db in DBS:
                    h = tax[db].get(mid)
                    if h:
                        row.append(h["id"])
                        if db == "ncbi":
                            row.append(acc2taxid.get(h["id"], ""))
                        row += [h["pct"], h["genus"], h["species"],
                                h["lineage"]]
                    else:
                        row += [""] * (6 if db == "ncbi" else 5)
                w.writerow(row)
                n_rows += 1

    print(f"[taxonomy_table] {n_rows} microbes from {len(batches)} plates "
          f"-> {args.out}")


if __name__ == "__main__":
    main()
