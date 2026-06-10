#!/usr/bin/env bash
# isoTAX replica - global-alignment taxonomy vs NCBI 16S and SILVA.
#   vsearch --usearch_global  (Needleman-Wunsch, --iddef 2 = exclude terminal
#   gaps, isolateR's identity definition), top hit per query.
#
#   02_isotax.sh <query.fasta> <out_dir> <batch> [ref_dir]
set -euo pipefail

QUERY="$1"; OUT="$2"; BATCH="$3"; REF="${4:-data/ref}"
THREADS="${THREADS:-4}"
mkdir -p "$OUT"
HERE="$(cd "$(dirname "$0")" && pwd)"

run_vsearch () { # <db_fasta> <hits_tsv>
    vsearch --usearch_global "$QUERY" --db "$1" \
        --id 0.50 --iddef 2 --maxaccepts 1 --maxrejects 32 --top_hits_only \
        --threads "$THREADS" --strand both \
        --userfields query+target+id --userout "$2" --quiet
}

# --- NCBI 16S type-strain ----------------------------------------------------
if [[ -s "$REF/ncbi_16S.fasta" ]]; then
    run_vsearch "$REF/ncbi_16S.fasta" "$OUT/$BATCH.ncbi.hits.tsv"
    python "$HERE/02_assign.py" --hits "$OUT/$BATCH.ncbi.hits.tsv" \
        --db ncbi --acc2taxid "$REF/ncbi_16S.acc2taxid.tsv" \
        --out "$OUT/$BATCH.tax.ncbi.csv"
else
    echo "  skip NCBI: $REF/ncbi_16S.fasta missing (run setup_db.sh)" >&2
fi

# --- SILVA SSU NR99 ----------------------------------------------------------
if [[ -s "$REF/silva_ssu.fasta" ]]; then
    run_vsearch "$REF/silva_ssu.fasta" "$OUT/$BATCH.silva.hits.tsv"
    python "$HERE/02_assign.py" --hits "$OUT/$BATCH.silva.hits.tsv" \
        --db silva --silva "$REF/silva_ssu.fasta" \
        --out "$OUT/$BATCH.tax.silva.csv"
else
    echo "  skip SILVA: $REF/silva_ssu.fasta missing (run setup_db.sh)" >&2
fi
