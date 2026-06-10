#!/usr/bin/env python3
"""Join SILVA and GTDB per-read taxonomy and score name concordance.

GTDB appends polyphyly suffixes to genus/species (e.g. Bacillus_A,
Pseudomonas_E). A "loose" comparison strips a trailing _<UPPER> so that
GTDB 'Bacillus_A' matches SILVA 'Bacillus'.

Usage:
  05_join.py --silva <tax.silva.csv> --gtdb <tax.gtdb.csv> \
             --out <combined.csv> --concordance <concordance.csv>
Prints a one-line summary; exits 0 even if one side is empty.
"""
import argparse
import csv
import re
import sys

SUFFIX = re.compile(r"_[A-Z]+$")          # GTDB polyphyly suffix, e.g. _A, _BA
# SILVA placeholder "species" that are not real binomials
PLACEHOLDER = {"sp.", "sp", "cf.", "aff.", "uncultured", "unidentified",
               "bacterium", "metagenome", "endophyte", "spp."}


def strip_suffix(tok):
    return SUFFIX.sub("", tok)


def loose_genus(name):
    """Genus with GTDB polyphyly suffix removed; placeholders -> ''."""
    g = name.strip()
    if not g or g.lower() in PLACEHOLDER:
        return ""
    return strip_suffix(g)


def species_epithet(name):
    """Lowercase species epithet, GTDB suffix removed; '' if not a real species.
    Handles 'Bacillus_A cereus_BA' -> 'cereus', 'Pseudomonas sp. TK30' -> ''."""
    if not name:
        return ""
    toks = name.split()
    if len(toks) < 2:
        return ""
    if toks[0].lower() in PLACEHOLDER:        # e.g. 'uncultured bacterium'
        return ""
    epi = strip_suffix(toks[1]).lower()
    return "" if epi in PLACEHOLDER else epi


def read_tax(path):
    """query -> (pct_id, genus, species)."""
    out = {}
    try:
        with open(path) as fh:
            r = csv.DictReader(fh)
            for row in r:
                out[row["query"]] = (row.get("pct_id", ""),
                                     row.get("genus", "").strip(),
                                     row.get("species", "").strip())
    except FileNotFoundError:
        pass
    return out


def cat(a, b):
    """Concordance category for two names (already comparable)."""
    if not a and not b:
        return "both_empty"
    if not a or not b:
        return "one_empty"
    return "match" if a == b else "mismatch"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--silva", required=True)
    ap.add_argument("--gtdb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--concordance", required=True)
    args = ap.parse_args()

    silva = read_tax(args.silva)
    gtdb = read_tax(args.gtdb)
    queries = sorted(set(silva) | set(gtdb))

    tally = {"genus_strict": {}, "genus_loose": {}, "species_epithet": {}}

    def bump(key, c):
        tally[key][c] = tally[key].get(c, 0) + 1

    with open(args.out, "w", newline="") as fo, \
            open(args.concordance, "w", newline="") as fc:
        wo = csv.writer(fo)
        wc = csv.writer(fc)
        wo.writerow(["query",
                     "silva_pct", "silva_genus", "silva_species",
                     "gtdb_pct", "gtdb_genus", "gtdb_species"])
        wc.writerow(["query", "silva_genus", "gtdb_genus",
                     "genus_strict", "genus_loose",
                     "silva_epithet", "gtdb_epithet", "species_epithet"])
        for q in queries:
            sp, sg, ss = silva.get(q, ("", "", ""))
            gp, gg, gs = gtdb.get(q, ("", "", ""))
            wo.writerow([q, sp, sg, ss, gp, gg, gs])

            gstrict = cat(sg, gg)
            gloose = cat(loose_genus(sg), loose_genus(gg))
            se, ge = species_epithet(ss), species_epithet(gs)
            sepi = cat(se, ge)
            bump("genus_strict", gstrict)
            bump("genus_loose", gloose)
            bump("species_epithet", sepi)
            wc.writerow([q, sg, gg, gstrict, gloose, se, ge, sepi])

    def pct(key, c):
        t = sum(tally[key].values()) or 1
        return f"{tally[key].get(c,0)}/{t}"

    print(f"[silva_vs_gtdb] {len(queries)} reads  "
          f"genus match(loose)={pct('genus_loose','match')}  "
          f"genus mismatch(loose)={pct('genus_loose','mismatch')}  "
          f"species match(epithet)={pct('species_epithet','match')}  "
          f"-> {args.out}")


if __name__ == "__main__":
    main()
