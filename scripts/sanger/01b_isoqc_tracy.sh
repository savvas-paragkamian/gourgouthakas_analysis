#!/usr/bin/env bash
# isoQC via tracy (pure-CLI alternative to 01_isoqc.py) + a comparison table.
#
#   01b_isoqc_tracy.sh <ab1_dir> <out_dir> <batch>
#
# Produces:
#   <out_dir>/<batch>.tracy.fasta   trimmed reads (tracy basecall)
#   <out_dir>/<batch>.tracy.csv     id,tracy_len
#   <out_dir>/<batch>.compare.csv   id,biopython_len,tracy_len  (if QC csv exists)
set -euo pipefail

ABDIR="$1"; OUT="$2"; BATCH="$3"
mkdir -p "$OUT" "$OUT/.tracy_tmp"
FA="$OUT/$BATCH.tracy.fasta"
TCSV="$OUT/$BATCH.tracy.csv"
: > "$FA"
echo "id,tracy_len" > "$TCSV"

shopt -s nullglob
for ab1 in "$ABDIR"/*.ab1; do
    id=$(basename "$ab1" .ab1)
    pref="$OUT/.tracy_tmp/$id"
    # tracy writes to exactly the -o path (it does NOT append .fasta).
    tracy basecall -f fasta -y primary -o "$pref.fasta" "$ab1" >/dev/null 2>&1 || { echo "  tracy failed: $id" >&2; continue; }
    if [[ -s "$pref.fasta" ]]; then
        # normalise header to the read id, append to batch fasta
        seqkit replace -p '.*' -r "$id" "$pref.fasta" >> "$FA" 2>/dev/null
        len=$(seqkit fx2tab -nl "$pref.fasta" 2>/dev/null | awk '{print $2; exit}')
        echo "$id,${len:-0}" >> "$TCSV"
    fi
done
rm -rf "$OUT/.tracy_tmp"
echo "[isoQC:tracy] -> $FA"

# Compare against the Biopython QC table if present.
QC="$OUT/$BATCH.qc.csv"
if [[ -f "$QC" ]]; then
    CMP="$OUT/$BATCH.compare.csv"
    awk -F, 'NR>1{bp[$1]=$3} END{}' "$QC" > /dev/null
    {
        echo "id,biopython_len,tracy_len"
        awk -F, 'FNR==NR{ if(FNR>1) bp[$1]=$3; next }
                 FNR>1{ print $1","(($1 in bp)?bp[$1]:"NA")","$2 }' \
            "$QC" "$TCSV"
    } > "$CMP"
    echo "[compare] -> $CMP"
fi
