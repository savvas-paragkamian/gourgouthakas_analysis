#!/usr/bin/env python3
"""Merge the SILVA, NCBI and GTDB taxonomy of every isolate into one TSV, join
the cave-isolate metadata, and build the per-depth abundance table.

Reads the per-batch isoTAX outputs (results/sanger/<batch>/<batch>.tax.<db>.csv)
rather than the merged ones: read ids restart per plate (e.g. 1_27f-A01 exists
on several plates) and the merged step keeps only the first of each collision,
so the per-batch files are the complete, unambiguous source. Each microbe is
keyed by (plate, microbe_id); for every database the reference identifier
(best hit), % identity and lineage are reported side by side.

The numeric "stab" id of each isolate is the leading number of microbe_id
(1280_27f-F12 -> 1280, 1056-Premixed -> 1056). It is inner-joined against
data/gourgouthakas-cave-isolates.csv (tab-separated): only isolates present in
that sheet are kept, each carrying its sampling metadata -- crucially the
`depth` column.

From that joined table the per-depth abundance matrices used by the tree figures
(06_gtdb_tree.R / 07_fasttree_genus_tree.R) are written, one per taxonomy:
results/gourgouthakas_depth_table.gtdb.tsv and .silva.tsv. Each has one row per
taxon (that taxonomy's genus and species labels) and one column per Gourgouthakas
sampling depth, counting the isolates of that taxon found at that depth. The
tables are kept separate so a figure only ever sees its own taxonomy -- mixing
them would, e.g., count the GTDB Aquipseudomonas isolates again under "Pseudomonas"
(their SILVA genus) in the GTDB figure. The CSV `depth` holds the positive
magnitude (m); the table columns are the signed depths (0, -39, ... -1100).

Usage:
  08_taxonomy_table.py [--results-dir results/sanger] \
                       [--out results/taxonomy_per_microbe.tsv] \
                       [--isolates data/gourgouthakas-cave-isolates.csv] \
                       [--depth-out results/gourgouthakas_depth_table.tsv]
"""
import argparse
import csv
import glob
import os
import re
from collections import defaultdict

DBS = ["silva", "ncbi", "gtdb"]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]

# Databases whose genus/species names are used as taxon labels in the tree
# figures (and therefore as rows of the depth table). NCBI is not plotted.
DEPTH_LABEL_DBS = ["gtdb", "silva"]

# The 9 Gourgouthakas sampling depths (metres). The metadata CSV stores the
# positive magnitude; the figures expect the signed value as the column name.
DEPTH_MAGNITUDES = [0, 39, 220, 418, 678, 713, 900, 1050, 1100]
DEPTH_COLS = ["0"] + [f"-{m}" for m in DEPTH_MAGNITUDES[1:]]
MAG_TO_COL = {0: "0", **{m: f"-{m}" for m in DEPTH_MAGNITUDES[1:]}}


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


def stab_of(microbe_id):
    """Leading numeric id shared with the metadata sheet (None if absent)."""
    m = re.match(r"(\d+)", microbe_id)
    return m.group(1) if m else None


def read_isolates(path):
    """stab -> metadata dict; plus the ordered list of metadata columns."""
    if not os.path.isfile(path):
        return {}, []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = [c for c in (reader.fieldnames or []) if c != "stab"]
        meta = {}
        for row in reader:
            key = (row.get("stab") or "").strip()
            if key:
                meta[key] = {c: (row.get(c) or "").strip() for c in cols}
    return meta, cols


def depth_col(value):
    """Map a CSV depth magnitude to its signed column name, or None."""
    try:
        mag = int(float(value))
    except (TypeError, ValueError):
        return None
    return MAG_TO_COL.get(mag)


def depth_path(base, db):
    """Per-taxonomy depth-table path: <base>.tsv -> <base>.<db>.tsv."""
    root, ext = os.path.splitext(base)
    return f"{root}.{db}{ext}"


def write_depth_table(counts, path):
    """taxon x depth count matrix -> TSV (backing up any existing file)."""
    if not counts:
        print(f"[taxonomy_table] no isolates with a known depth; "
              f"left {path} untouched")
        return
    backed_up = os.path.isfile(path)
    if backed_up:
        os.replace(path, path + ".bak")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["taxon"] + DEPTH_COLS)
        for taxon in sorted(counts):
            w.writerow([taxon] + [counts[taxon].get(c, 0) for c in DEPTH_COLS])
    print(f"[taxonomy_table] {len(counts)} taxa x {len(DEPTH_COLS)} depths "
          f"-> {path}" + (f" (previous -> {os.path.basename(path)}.bak)"
                          if backed_up else ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results/sanger")
    ap.add_argument("--out", default="results/taxonomy_per_microbe.tsv")
    ap.add_argument("--acc2taxid", default="data/ref/ncbi_16S.acc2taxid.tsv",
                    help="NCBI accession<TAB>taxid map (from setup_db.sh)")
    ap.add_argument("--isolates",
                    default="data/gourgouthakas-cave-isolates.csv",
                    help="tab-separated isolate metadata keyed by `stab`")
    ap.add_argument("--depth-out",
                    default="results/gourgouthakas_depth_table.tsv",
                    help="base path for the per-taxonomy depth tables; the "
                         "suffix .<db>.tsv is inserted (.gtdb.tsv / .silva.tsv)")
    args = ap.parse_args()

    # NCBI accession -> taxid (the best_hit in *.tax.ncbi.csv is the accession)
    acc2taxid = {}
    if os.path.isfile(args.acc2taxid):
        with open(args.acc2taxid) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    acc2taxid[p[0]] = p[1]

    isolates, meta_cols = read_isolates(args.isolates)
    if isolates:
        print(f"[taxonomy_table] {len(isolates)} isolates with metadata from "
              f"{args.isolates}")
    else:
        print(f"[taxonomy_table] no metadata loaded from {args.isolates}; "
              f"depth table will be empty")

    batches = sorted(
        d for d in glob.glob(os.path.join(args.results_dir, "*"))
        if os.path.isdir(d) and os.path.basename(d) != "merged")

    header = ["plate", "microbe_id", "stab"]
    for db in DBS:
        header.append(f"{db}_id")
        if db == "ncbi":
            header.append("ncbi_taxid")
        header += [f"{db}_pct_id", f"{db}_genus", f"{db}_species",
                   f"{db}_lineage"]
    header += meta_cols

    # one depth table per taxonomy (db -> taxon -> depth column -> count).
    # Keeping them separate avoids mixing labels: e.g. the GTDB split of the
    # Pseudomonas complex (Aquipseudomonas, Stutzerimonas, ...) must not pick up
    # the isolates SILVA still calls "Pseudomonas", which would double-count them
    # against the GTDB genera in the GTDB figure.
    depth_counts = {db: defaultdict(lambda: defaultdict(int))
                    for db in DEPTH_LABEL_DBS}
    n_with_depth = 0

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
                stab = stab_of(mid)
                meta = isolates.get(stab) if stab else None
                if meta is None:
                    continue   # inner join: only isolates in the metadata sheet
                row = [plate, mid, stab or ""]
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
                row += [meta.get(c, "") for c in meta_cols]
                w.writerow(row)
                n_rows += 1

                # per-depth taxon counts (every retained isolate has metadata),
                # accumulated separately per taxonomy
                col = depth_col(meta.get("depth", ""))
                if col is not None:
                    n_with_depth += 1
                    for db in DEPTH_LABEL_DBS:
                        h = tax[db].get(mid)
                        if not h:
                            continue
                        for lab in {h["genus"], h["species"]}:
                            if lab:
                                depth_counts[db][lab][col] += 1

    print(f"[taxonomy_table] {n_rows} metadata-matched microbes (inner join) "
          f"from {len(batches)} plates -> {args.out}; "
          f"{n_with_depth} with a usable depth")

    for db in DEPTH_LABEL_DBS:
        write_depth_table(depth_counts[db], depth_path(args.depth_out, db))


if __name__ == "__main__":
    main()
