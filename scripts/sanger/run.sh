#!/usr/bin/env bash
# isolateR-equivalent Sanger 16S pipeline (isoQC -> isoTAX -> isoLIB -> tree).
# Runs inside the sanger16s container. Mount the repo at /work.
#
#   podman run --rm -v "$PWD":/work -w /work sanger16s bash scripts/sanger/run.sh
#
# Env overrides:
#   QC_MODE=auto|window   (default auto; isolateR "auto cutoff")
#   THREADS=4
#   REF=data/ref          reference DB dir
#   OUT=results/sanger     output dir
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REF="${REF:-data/ref}"
OUT="${OUT:-results/sanger}"
QC_MODE="${QC_MODE:-auto}"
export THREADS="${THREADS:-4}"
mkdir -p "$OUT"

# Input batches = every ABI plate dir + the Villy premixed set.
mapfile -t BATCHDIRS < <(find data -maxdepth 1 -type d \
    \( -name '*ABI' -o -name '*_ab1' \) | sort)

if [[ ${#BATCHDIRS[@]} -eq 0 ]]; then
    echo "No ABI directories found under data/" >&2; exit 1
fi

# Ensure reference DBs exist (downloads once).
if [[ ! -s "$REF/ncbi_16S.fasta" || ! -s "$REF/silva_ssu.fasta" ]]; then
    echo "== reference DBs not found, running setup_db.sh =="
    bash "$HERE/setup_db.sh" "$REF"
fi
export TAXONKIT_DB="$PWD/$REF/taxdump"

ALL_QC="$OUT/all.qc.fasta"; : > "$ALL_QC"

for dir in "${BATCHDIRS[@]}"; do
    batch="$(basename "$dir")"
    bout="$OUT/$batch"; mkdir -p "$bout"
    echo "==================== $batch ===================="

    # 1. isoQC (Biopython, auto/Mott by default) + tracy comparison
    python "$HERE/01_isoqc.py" --in "$dir" --mode "$QC_MODE" \
        --fasta "$bout/$batch.qc.fasta" --csv "$bout/$batch.qc.csv"
    bash "$HERE/01b_isoqc_tracy.sh" "$dir" "$bout" "$batch" || true

    [[ -s "$bout/$batch.qc.fasta" ]] || { echo "  no PASS reads, skipping"; continue; }
    cat "$bout/$batch.qc.fasta" >> "$ALL_QC"

    # 2. isoTAX (NCBI + SILVA)
    bash "$HERE/02_isotax.sh" "$bout/$batch.qc.fasta" "$bout" "$batch" "$REF"
    # 3. isoLIB
    bash "$HERE/03_isolib.sh" "$bout/$batch.qc.fasta" "$bout" "$batch"
    # 4. tree
    bash "$HERE/04_tree.sh" "$bout/$batch.reps.fasta" "$bout" "$batch"
done

# ---- merged cross-plate library --------------------------------------------
# Same order/input as per-batch: classify the reads first (plain ids, no
# vsearch ;size= suffix) so 03_isolib.sh can join taxonomy onto each rep.
echo "==================== merged ===================="
mout="$OUT/merged"; mkdir -p "$mout"
cp "$ALL_QC" "$mout/merged.qc.fasta"
bash "$HERE/02_isotax.sh"  "$mout/merged.qc.fasta"  "$mout" merged "$REF"
bash "$HERE/03_isolib.sh"  "$mout/merged.qc.fasta"  "$mout" merged
bash "$HERE/04_tree.sh"    "$mout/merged.reps.fasta" "$mout" merged

echo "== done. results under $OUT/ =="
