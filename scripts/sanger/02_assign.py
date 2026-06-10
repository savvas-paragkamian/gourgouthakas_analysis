#!/usr/bin/env python3
"""Attach lineage to vsearch hits and assign the finest confident rank.

isolateR isoTAX rank cutoffs (16S, % identity):
  phylum 75.0  class 78.5  order 82.0  family 86.5  genus 96.5  species 98.7
The assigned rank is the finest rank whose cutoff the hit identity meets;
the reported lineage is truncated to that rank.

Usage:
  02_assign.py --hits hits.tsv --db ncbi  --acc2taxid map.tsv  --out tax.csv
  02_assign.py --hits hits.tsv --db silva --silva silva.fasta  --out tax.csv
  02_assign.py --hits hits.tsv --db gtdb  --gtdb gtdb_ssu.fasta --out tax.csv

hits.tsv columns (vsearch --userfields query+target+id): query, target, id
"""
import argparse
import csv
import subprocess
import sys

RANKS = ["phylum", "class", "order", "family", "genus", "species"]
CUTOFF = {"phylum": 75.0, "class": 78.5, "order": 82.0,
          "family": 86.5, "genus": 96.5, "species": 98.7}


def read_hits(path):
    hits = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[0] not in hits:  # keep top hit only
                hits[p[0]] = (p[1], float(p[2]))
    return hits


def lineages_ncbi(targets, acc2taxid):
    """target accession -> [phylum..species] via taxonkit."""
    a2t = {}
    with open(acc2taxid) as fh:
        for line in fh:
            a, t = line.rstrip("\n").split("\t")[:2]
            a2t[a] = t
    taxids = sorted({a2t[t] for t in targets if t in a2t})
    if not taxids:
        return {}, a2t
    # taxonkit lineage -R gives lineage + matching ranks
    proc = subprocess.run(
        ["taxonkit", "lineage", "-R"],
        input="\n".join(taxids) + "\n",
        text=True, capture_output=True, check=True)
    by_taxid = {}
    for line in proc.stdout.splitlines():
        c = line.split("\t")
        if len(c) < 3:
            continue
        taxid, names, ranks = c[0], c[1].split(";"), c[2].split(";")
        rmap = {r: n for n, r in zip(names, ranks)}
        by_taxid[taxid] = [rmap.get(r, "") for r in RANKS]
    return {t: by_taxid.get(a2t.get(t, ""), [""] * len(RANKS))
            for t in targets}, a2t


def lineages_silva(targets, silva_fasta):
    """SILVA seqid -> [phylum..species]. Header: '>ID Domain;Phylum;...;Genus species'."""
    want = set(targets)
    found = {}
    with open(silva_fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            head = line[1:].rstrip("\n")
            sid = head.split()[0]
            if sid not in want:
                continue
            tax = head.split(" ", 1)[1] if " " in head else ""
            parts = [x.strip() for x in tax.split(";") if x.strip()]
            # SILVA path: Domain;Phylum;Class;Order;Family;Genus(;species)
            body = parts[1:] if parts else []          # drop Domain
            lin = (body + [""] * len(RANKS))[:len(RANKS)]
            found[sid] = lin
            if len(found) == len(want):
                break
    return {t: found.get(t, [""] * len(RANKS)) for t in targets}


def lineages_gtdb(targets, gtdb_fasta):
    """GTDB seqid -> [phylum..species]. Header:
    '>GCA_x d__..;p__..;c__..;o__..;f__..;g__..;s__Name sp [locus_tag=..]'."""
    want = set(targets)
    pref = {"p": "phylum", "c": "class", "o": "order",
            "f": "family", "g": "genus", "s": "species"}
    found = {}
    with open(gtdb_fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            head = line[1:].rstrip("\n")
            sid = head.split(" ", 1)[0]
            if sid not in want:
                continue
            rest = head.split(" ", 1)[1] if " " in head else ""
            tax = rest.split(" [", 1)[0]          # strip trailing [locus_tag=..]
            rmap = {}
            for tok in tax.split(";"):
                tok = tok.strip()
                if len(tok) > 3 and tok[1:3] == "__" and tok[0] in pref:
                    rmap[pref[tok[0]]] = tok[3:]
            found[sid] = [rmap.get(r, "") for r in RANKS]
            if len(found) == len(want):
                break
    return {t: found.get(t, [""] * len(RANKS)) for t in targets}


def assign(lin, ident):
    """Truncate lineage to the finest rank the identity supports."""
    best = ""
    out = [""] * len(RANKS)
    for i, r in enumerate(RANKS):
        if ident >= CUTOFF[r] and lin[i]:
            out[i] = lin[i]
            best = r
    return best, ";".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True)
    ap.add_argument("--db", choices=["ncbi", "silva", "gtdb"], required=True)
    ap.add_argument("--acc2taxid")
    ap.add_argument("--silva")
    ap.add_argument("--gtdb")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    hits = read_hits(args.hits)
    targets = [t for t, _ in hits.values()]

    if args.db == "ncbi":
        if not args.acc2taxid:
            sys.exit("--acc2taxid required for --db ncbi")
        lin_by_target, _ = lineages_ncbi(targets, args.acc2taxid)
    elif args.db == "silva":
        if not args.silva:
            sys.exit("--silva required for --db silva")
        lin_by_target = lineages_silva(targets, args.silva)
    else:
        if not args.gtdb:
            sys.exit("--gtdb required for --db gtdb")
        lin_by_target = lineages_gtdb(targets, args.gtdb)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["query", "best_hit", "pct_id", "assigned_rank",
                    "phylum", "class", "order", "family", "genus", "species"])
        for q in sorted(hits):
            target, ident = hits[q]
            lin = lin_by_target.get(target, [""] * len(RANKS))
            rank, _ = assign(lin, ident)
            shown = [lin[i] if ident >= CUTOFF[RANKS[i]] else ""
                     for i in range(len(RANKS))]
            w.writerow([q, target, ident, rank, *shown])
    print(f"[isoTAX:{args.db}] {len(hits)} assignments -> {args.out}")


if __name__ == "__main__":
    main()
