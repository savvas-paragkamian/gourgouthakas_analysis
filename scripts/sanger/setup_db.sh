#!/usr/bin/env bash
# Download reference databases into data/ref/ (run once; idempotent).
#   - NCBI 16S_ribosomal_RNA  == RefSeq Targeted Loci type material
#                                (isolateR's default isoTAX database)
#   - NCBI taxdump            == lineages for taxonkit
#   - SILVA SSU Ref NR99      == broader 16S/18S reference
#
#   setup_db.sh [ref_dir]      (default: data/ref)
set -euo pipefail

REF="${1:-data/ref}"
mkdir -p "$REF"
cd "$REF"

SILVA_URL="https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"

# ---- NCBI 16S type-strain DB -> FASTA (with taxids) -------------------------
# Fetch the prebuilt BLAST DB tarball directly with wget (update_blastdb.pl
# requires curl, which is not in the image).
if [[ ! -s ncbi_16S.fasta ]]; then
    echo "[db] downloading NCBI 16S_ribosomal_RNA ..."
    mkdir -p blastdb && cd blastdb
    wget -q https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
    wget -q https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
    tar -xzf 16S_ribosomal_RNA.tar.gz && tar -xzf taxdb.tar.gz
    rm -f 16S_ribosomal_RNA.tar.gz taxdb.tar.gz
    cd ..
    echo "[db] exporting FASTA + acc->taxid map ..."
    blastdbcmd -db blastdb/16S_ribosomal_RNA -entry all -outfmt "%f" > ncbi_16S.fasta
    blastdbcmd -db blastdb/16S_ribosomal_RNA -entry all -outfmt "%a	%T" > ncbi_16S.acc2taxid.tsv
fi

# ---- NCBI taxdump for taxonkit ---------------------------------------------
if [[ ! -s taxdump/nodes.dmp ]]; then
    echo "[db] downloading NCBI taxdump ..."
    wget -q https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    mkdir -p taxdump
    tar -xzf taxdump.tar.gz -C taxdump nodes.dmp names.dmp delnodes.dmp merged.dmp
    rm -f taxdump.tar.gz
fi

# ---- SILVA SSU Ref NR99 -----------------------------------------------------
if [[ ! -s silva_ssu.fasta ]]; then
    echo "[db] downloading SILVA SSU Ref NR99 (138.2) ..."
    wget -q -O silva_ssu.fasta.gz "$SILVA_URL"
    gzip -dc silva_ssu.fasta.gz | sed '/^>/!s/U/T/g; /^>/!s/u/t/g' > silva_ssu.fasta
    rm -f silva_ssu.fasta.gz
fi

# ---- GTDB SSU representative sequences --------------------------------------
# Bacterial + archaeal SSU from GTDB species-representative genomes, carrying
# GTDB taxonomy in the header (d__;p__;c__;o__;f__;g__;s__).
GTDB_REL="https://data.gtdb.ecogenomic.org/releases/latest"
GTDB_BASE="$GTDB_REL/genomic_files_reps"
if [[ ! -s gtdb_ssu.fasta ]]; then
    echo "[db] downloading GTDB SSU reps (bac120 + ar53) ..."
    wget -q -O bac120_ssu_reps.fna.gz "$GTDB_BASE/bac120_ssu_reps.fna.gz"
    wget -q -O ar53_ssu_reps.fna.gz   "$GTDB_BASE/ar53_ssu_reps.fna.gz"
    gzip -dc bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz > gtdb_ssu.fasta
    rm -f bac120_ssu_reps.fna.gz ar53_ssu_reps.fna.gz
fi

# ---- GTDB master tree + taxonomy (for pruning in 06_gtdb_tree.R) ------------
mkdir -p gtdb
if [[ ! -s gtdb/bac120.tree ]]; then
    echo "[db] downloading GTDB bac120 master tree ..."
    wget -q -O gtdb/bac120.tree.gz "$GTDB_REL/bac120.tree.gz"
    gzip -df gtdb/bac120.tree.gz
fi
if [[ ! -s gtdb/bac120_taxonomy.tsv ]]; then
    echo "[db] downloading GTDB bac120 taxonomy ..."
    wget -q -O gtdb/bac120_taxonomy.tsv.gz "$GTDB_REL/bac120_taxonomy.tsv.gz"
    gzip -df gtdb/bac120_taxonomy.tsv.gz
fi

echo "[db] ready:"
ls -lh ncbi_16S.fasta silva_ssu.fasta gtdb_ssu.fasta \
       gtdb/bac120.tree gtdb/bac120_taxonomy.tsv taxdump/nodes.dmp 2>/dev/null
