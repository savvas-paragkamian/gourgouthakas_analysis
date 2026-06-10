#!/usr/bin/env bash
# isoLIB replica - dereplicate PASS reads into a strain library at 99.5% id.
#   vsearch --cluster_size --id 0.995  (isolateR default group cutoff)
#
#   03_isolib.sh <query.fasta> <out_dir> <batch>
# Outputs:
#   <batch>.reps.fasta   cluster representatives (centroids)
#   <batch>.clusters.uc  raw vsearch clustering
#   <batch>.lib.csv      rep, cluster_size, members, rep taxonomy (NCBI if present)
set -euo pipefail

QUERY="$1"; OUT="$2"; BATCH="$3"
THREADS="${THREADS:-4}"
mkdir -p "$OUT"

vsearch --cluster_size "$QUERY" --id 0.995 --iddef 2 --strand both \
    --threads "$THREADS" --sizeout \
    --centroids "$OUT/$BATCH.reps.fasta" \
    --uc "$OUT/$BATCH.clusters.uc" --quiet

TAX="$OUT/$BATCH.tax.ncbi.csv"
LIB="$OUT/$BATCH.lib.csv"

# members per centroid from the .uc (H = member, S = seed/centroid)
awk -F'\t' -v tax="$TAX" '
    BEGIN{
        if (tax != "" && (getline _ < tax) > 0) {              # header line
            while ((getline l < tax) > 0) {
                n=split(l,a,","); tx[a[1]]=a[4]"|"a[9]" "a[10] # rank|genus species
            }
        }
    }
    $1=="S"{ rep[$9]=$9; size[$9]=1; members[$9]=$9 }
    $1=="H"{ size[$10]++; members[$10]=members[$10]";"$9 }
    END{
        print "representative,cluster_size,rep_taxonomy,members"
        for (r in size)
            printf "%s,%d,%s,%s\n", r, size[r], (r in tx?tx[r]:"NA"), members[r]
    }' "$OUT/$BATCH.clusters.uc" | sort -t, -k2,2nr > "$LIB"

NREP=$(grep -c "^>" "$OUT/$BATCH.reps.fasta" || echo 0)
NQ=$(grep -c "^>" "$QUERY" || echo 0)
echo "[isoLIB] $NQ reads -> $NREP strain clusters -> $LIB"
