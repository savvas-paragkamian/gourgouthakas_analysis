#!/usr/bin/env bash
# Phylogeny of strain-library representatives.
#   mafft --auto  ->  FastTree -nt -gtr
#
#   04_tree.sh <reps.fasta> <out_dir> <batch>
set -euo pipefail

REPS="$1"; OUT="$2"; BATCH="$3"
mkdir -p "$OUT"

n=$(grep -c "^>" "$REPS" || echo 0)
if [[ "$n" -lt 3 ]]; then
    echo "[tree] only $n reps; need >=3 for a tree, skipping" >&2
    exit 0
fi

# Strip the vsearch ';size=N' suffix from labels first: ';' is the Newick
# terminator, so leaving it in produces an unparseable tree.
CLEAN="$OUT/$BATCH.reps.clean.fasta"
seqkit replace -p ';size=[0-9]+' -r '' "$REPS" > "$CLEAN" 2>/dev/null

mafft --auto --thread "${THREADS:-4}" "$CLEAN" > "$OUT/$BATCH.reps.aln" 2>/dev/null
FastTree -nt -gtr "$OUT/$BATCH.reps.aln" > "$OUT/$BATCH.reps.nwk" 2>/dev/null
rm -f "$CLEAN"
echo "[tree] $n reps -> $OUT/$BATCH.reps.nwk"
