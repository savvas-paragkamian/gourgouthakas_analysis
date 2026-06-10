#!/usr/bin/env bash
# Connect SILVA taxonomy with GTDB taxonomy.
#
# For every batch that already has SILVA results (results/sanger/<batch>/), this:
#   1. classifies the same reads against the GTDB SSU reps DB (vsearch, the same
#      --usearch_global approach as isoTAX)  ->  <batch>.tax.gtdb.csv
#   2. joins SILVA + GTDB per read                 ->  <batch>.silva_vs_gtdb.csv
#   3. scores genus/species name concordance       ->  <batch>.concordance.csv
# and writes a combined concordance summary across all batches.
#
#   bash scripts/sanger/05_silva_gtdb.sh            (run inside the container)
#
# Run AFTER run.sh. Standalone: it only needs the per-batch *.qc.fasta and
# *.tax.silva.csv, plus the GTDB DB (downloaded on demand into data/ref/).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REF="${REF:-data/ref}"
OUT="${OUT:-results/sanger}"
THREADS="${THREADS:-4}"

# Ensure the GTDB SSU DB exists (downloads bac120+ar53 reps once).
if [[ ! -s "$REF/gtdb_ssu.fasta" ]]; then
    echo "== GTDB SSU DB not found, fetching via setup_db.sh =="
    bash "$HERE/setup_db.sh" "$REF"
fi

SUMMARY="$OUT/silva_vs_gtdb.summary.csv"
echo "batch,reads,genus_match_loose,genus_mismatch_loose,species_match_epithet" > "$SUMMARY"

shopt -s nullglob
for silva in "$OUT"/*/*.tax.silva.csv; do
    bdir="$(dirname "$silva")"
    batch="$(basename "$bdir")"
    # locate this batch's query reads
    qc="$bdir/$batch.qc.fasta"
    [[ -f "$qc" ]] || qc="$bdir/merged.qc.fasta"          # merged batch
    [[ -f "$qc" ]] || { echo "  skip $batch: no qc.fasta" >&2; continue; }
    echo "==================== $batch ===================="

    # 1. classify reads vs GTDB
    vsearch --usearch_global "$qc" --db "$REF/gtdb_ssu.fasta" \
        --id 0.50 --iddef 2 --maxaccepts 1 --maxrejects 32 --top_hits_only \
        --threads "$THREADS" --strand both \
        --userfields query+target+id --userout "$bdir/$batch.gtdb.hits.tsv" --quiet
    python "$HERE/02_assign.py" --hits "$bdir/$batch.gtdb.hits.tsv" \
        --db gtdb --gtdb "$REF/gtdb_ssu.fasta" --out "$bdir/$batch.tax.gtdb.csv"

    # 2 + 3. join SILVA + GTDB and score concordance
    line="$(python "$HERE/05_join.py" \
        --silva "$silva" --gtdb "$bdir/$batch.tax.gtdb.csv" \
        --out "$bdir/$batch.silva_vs_gtdb.csv" \
        --concordance "$bdir/$batch.concordance.csv")"
    echo "$line"

    # pull the numbers (X/Y) out of the printed summary for the roll-up
    reads=$(awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/){print $2; exit}}' <<<"$line")
    gm=$(grep -oE 'genus match\(loose\)=[0-9]+/[0-9]+' <<<"$line" | cut -d= -f2)
    gx=$(grep -oE 'genus mismatch\(loose\)=[0-9]+/[0-9]+' <<<"$line" | cut -d= -f2)
    sm=$(grep -oE 'species match\(epithet\)=[0-9]+/[0-9]+' <<<"$line" | cut -d= -f2)
    echo "$batch,$reads,$gm,$gx,$sm" >> "$SUMMARY"
done

echo "== done. summary -> $SUMMARY =="
