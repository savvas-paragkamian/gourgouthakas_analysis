# Gourgouthakas cave analysis

Here we analyse the data from the biobank of microbial isolates from 
the Gourgouthakas cave. 

## Pangenome analysis with anvi'o

### 1) Pangenome analysis with publicly available genomes of *Streptomyces*

#### Get the publicly available genomes of selected *Streptomyces*

The publicly available genomes of selected *Streptomyces* were downloaded using the [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) command (v0.3.3). 

##### Make the working directory and the directory where the genomes will be downloaded 

```
mkdir Streptomyces_genomes_anvio
mkdir Streptomyces_genomes_anvio/Streptomyces_genomes
```

##### Download the genomes

```
ncbi-genome-download bacteria --assembly-accessions GCF_000009765.2,GCF_000203835.1,GCF_003248355.1,GCF_003287935.1,GCF_003287915.1,GCF_003054555.1,GCF_048572725.1,GCF_049310085.1,GCF_049310105.1,GCF_049310065.1,GCF_003248335.1 --formats fasta --flat-output --output-folder ./Streptomyces_genomes_anvio/Streptomyces_genomes 
gunzip ./Streptomyces_genomes_anvio/Streptomyces_genomes/*.fna.gz
```

#### Pangenome analysis 

1) Rename the files, change the contig naming in the fasta files and prepare a database for every genome

```
cd ./Streptomyces_genomes_anvio
mkdir ./Streptomyces_genomes_simplified 
mkdir ./Streptomyces_genomes_db
cd ./Streptomyces_genomes
```

```
for file in *.fna *.fasta; do
  [ -e "$file" ] || continue

  if [[ "$file" == GCF_*_genomic.fna ]]; then
    newname=$(echo "$file" | sed -E 's/^(GCF_[0-9]+)\.([0-9]+).*\.fna$/\1_\2.fna/')
    mv "$file" "$newname"

  elif [[ "$file" == [0-9]*.unicycler_hybrid_assembly.fasta ]]; then
    newname=$(echo "$file" | sed -E 's/^([0-9]+)\.unicycler_hybrid_assembly\.fasta$/SRL\1.fasta/')
    mv "$file" "$newname"
  fi
done
```

```
for fasta in *.fasta *.fna; do
  [ -e "$fasta" ] || continue
  base=$(basename "$fasta")
  prefix="${base%.*}"

  # Reformat fasta and simplify contig names
  anvi-script-reformat-fasta "$fasta" \
    -o "../Streptomyces_genomes_simplified/${prefix}_simplified.fna" \
    --seq-type NT \
    --simplify-names \
    --report-file "../Streptomyces_genomes_simplified/${prefix}_rename-report.txt" \
    --prefix "$prefix"
  
  anvi-gen-contigs-database -f "../Streptomyces_genomes_simplified/${prefix}_simplified.fna" -o "../Streptomyces_genomes_db/${prefix}.db" --project-name "$prefix" -T 20 

done
```

2) Annotate the db files

```
anvi-setup-kegg-data
anvi-setup-ncbi-cogs
anvi-setup-pfams
anvi-setup-scg-taxonomy
anvi-setup-trna-taxonomy
anvi-setup-interacdome
```
```
cd ../Streptomyces_genomes_db 
```
```
for db in *.db; do
  echo "🔍 Annotating $db"
  echo "🔍 1. HMMS Annotation"
  anvi-run-hmms -c "$db" -I Bacteria_71 --also-scan-trnas -T 20
  echo "🔍 2. NCBI COG Annotation"
  anvi-run-ncbi-cogs -c "$db" -T 20
  echo "🔍 3. KEGG Kofam Annotation"
  anvi-run-kegg-kofams -c "$db" -T 20
  echo "🔍 4. Pfam Annotation"
  anvi-run-pfams -c "$db" -T 20
  echo "🔍 5. SCG Annotation"
  anvi-run-scg-taxonomy --contigs-db "$db" -T 20
done
```

3) Create a genomes-storage database file

```
echo -e "name\tcontigs_db_path" > Streptomyces_genomes.txt
for db in *.db; do
  name=$(basename "$db" .db)
  echo -e "${name}\t$(realpath "$db")" >> Streptomyces_genomes.txt
done
```
```
anvi-gen-genomes-storage -e Streptomyces_genomes.txt -o Streptomyces-GENOMES.db
```

4) Run pangenome analysis

```
mkdir -p ../Streptomyces_pangenome
```

```
anvi-pan-genome \
  -g Streptomyces-GENOMES.db \
  -n Streptomyces_Pangenome \
  -o ../Streptomyces_pangenome/Streptomyces_Pangenome-PAN.db \
  -T 20 \
  --minbit 0.5 \
  --mcl-inflation 10 \
  --use-ncbi-blast
```

#### Extraction of gene cluster presence/absence data from the pangenome analysis

```
cd ../Streptomyces_pangenome
```

```
sqlite3 -header -separator $'\t' Streptomyces_Pangenome-PAN.db "SELECT gene_cluster_id, genome_name FROM gene_clusters;" > Streptomyces_gene_cluster_output.txt
```

## 2) Pangenome analysis of our three isolates 

```
cd ../Streptomyces_genomes_db
```

1) Create a genomes-storage database file for only the genomes of our three isolates

```
echo -e "name\tcontigs_db_path" > biobank_Streptomyces_genomes.txt
for db in SRL*.db; do
  name=$(basename "$db" .db)
  echo -e "${name}\t$(realpath "$db")" >> biobank_Streptomyces_genomes.txt
done
```
```
anvi-gen-genomes-storage -e biobank_Streptomyces_genomes.txt -o BiobankStreptomyces-GENOMES.db
```

2) Run pangenome analysis

```
mkdir -p ../BiobankStreptomyces_Pangenome
```

```
anvi-pan-genome \
  -g BiobankStreptomyces-GENOMES.db \
  -n BiobankStreptomyces_Pangenome \
  -o ../BiobankStreptomyces_Pangenome/BiobankStreptomyces_Pangenome-PAN.db \
  -T 20 \
  --minbit 0.5 \
  --mcl-inflation 10 \
  --use-ncbi-blast
```

### Extraction of gene cluster presence/absence data from the pangenome analysis

```
cd ../BiobankStreptomyces_Pangenome
```

```
sqlite3 -header -separator $'\t' BiobankStreptomyces_Pangenome-PAN.db "SELECT gene_cluster_id, genome_name FROM gene_clusters;" > BiobankStreptomyces_gene_cluster_output.txt
```

### Design of the Pangenome Graphs

A directory for the pangenome graphs was created in the "Streptomyces_genomes_anvio" directory, using the following command:

```
mkdir Streptomyces_genomes_anvio_Graphs
```

The files Streptomyces_gene_cluster_output.txt and BiobankStreptomyces_gene_cluster_output.txt were copied inside the "Streptomyces_genomes_anvio_Graphs". Then, the graphs were created in R (version 4.5.1) using the [scripts/pangenome_graphs.R](scripts/pangenome_graphs.R) script. This script was executed inside the "Streptomyces_genomes_anvio_Graphs" directory. 

## Obtain functional info about our three cave isolates from the pangenome analysis output

A directory for the output files of this analysis was created in the "Streptomyces_genomes_anvio" directory, using the following command:

```
mkdir Streptomyces_functions
```
```
cd Streptomyces_functions
```

1) Export the gene cluster membership table 

```
anvi-export-gene-clusters \
  -p ../Streptomyces_pangenome/Streptomyces_Pangenome-PAN.db \
  -o Streptomyces_gene_clusters.tsv
```

2) Export functional info for each of the three isolates

```
anvi-export-functions \
  -c ../Streptomyces_genomes_db/SRL740.db \
  -o SRL740_functions.tsv
```
```
anvi-export-functions \
  -c ../Streptomyces_genomes_db/SRL742.db \
  -o SRL742_functions.tsv
```
```
anvi-export-functions \
  -c ../Streptomyces_genomes_db/SRL1060.db \
  -o SRL1060_functions.tsv
```

```
awk 'FNR==1 && NR!=1 {next} {print}' SRL1060_functions.tsv SRL740_functions.tsv SRL742_functions.tsv \
  > SRL_isolates_functions.tsv
```

 The exported tables were analyzed in R (version 4.5.1) using the [pangenome_functions.R](scripts/pangenome_functions.R) script. This script was executed inside the "Streptomyces_functions" directory.

## Taxonomic identification of isolates from Sanger 16S sequencing

The microbial isolates were identified by Sanger sequencing of the 16S rRNA gene
(27f primer). The raw `.ab1` trace files are organised under `data/` as one
directory per sequencing plate, plus the `Villy-Katerina_90-1068815155_ab1`
premixed set.

These were processed with a reproduction of the [isolateR](https://github.com/bdaisley/isolateR)
workflow (Daisley et al. 2024, *Bioinformatics* 40(7):btae448). Because the R
package could not be installed, the workflow was rebuilt with command-line tools
running inside a [podman](https://podman.io/) container, defined in
[scripts/sanger/Containerfile](scripts/sanger/Containerfile). The container
provides the bioconda toolchain (`vsearch`, `tracy`, `seqkit`, `blast`,
`taxonkit`, `mafft`, `fasttree`, `biopython`). The defaults mirror isolateR:
auto quality-trim cutoff (or fixed Phred 20 / window 15 / minimum length 200 bp),
identity excluding terminal gaps (`vsearch --iddef 2`), taxonomic rank cutoffs
(phylum 75 / class 78.5 / order 82 / family 86.5 / genus 96.5 / species 98.7 %),
and strain-library dereplication at 99.5 % identity.

| isolateR function | Step | Tool in this pipeline |
|---|---|---|
| `isoQC`  | quality-trim `.ab1` → FASTA + PASS/FAIL table | [01_isoqc.py](scripts/sanger/01_isoqc.py) (Biopython) and [01b_isoqc_tracy.sh](scripts/sanger/01b_isoqc_tracy.sh) (tracy), compared |
| `isoTAX` | taxonomy by global alignment | [02_isotax.sh](scripts/sanger/02_isotax.sh) → `vsearch --usearch_global` vs **NCBI 16S** type strains and **SILVA** SSU NR99, lineages via `taxonkit` ([02_assign.py](scripts/sanger/02_assign.py)) |
| `isoLIB` | strain library | [03_isolib.sh](scripts/sanger/03_isolib.sh) → `vsearch --cluster_size --id 0.995` |
| (optional) | phylogeny of representatives | [04_tree.sh](scripts/sanger/04_tree.sh) → `mafft` + `FastTree` |

Build the image and run the full pipeline from the repository root:

```
podman build -t sanger16s -f scripts/sanger/Containerfile .
podman run --rm -v "$PWD":/work -w /work sanger16s bash scripts/sanger/run.sh
```

The reference databases are downloaded once into `data/ref/` by
[scripts/sanger/setup_db.sh](scripts/sanger/setup_db.sh) (NCBI 16S RefSeq
Targeted Loci, NCBI taxdump, SILVA SSU Ref NR99 release 138.2). The orchestrator
[scripts/sanger/run.sh](scripts/sanger/run.sh) processes each plate as a separate
batch and then builds a merged, cross-plate strain library. Results are written
to `results/sanger/<batch>/`:

```
<batch>.qc.fasta / .qc.csv          isoQC trimmed reads + PASS/FAIL table
<batch>.tracy.fasta / .compare.csv  tracy trimming + length comparison
<batch>.tax.ncbi.csv                taxonomy vs NCBI 16S type strains
<batch>.tax.silva.csv               taxonomy vs SILVA SSU NR99
<batch>.reps.fasta / .lib.csv       strain library (99.5 % clusters)
<batch>.reps.aln / .reps.nwk        alignment and phylogenetic tree
```

The behaviour can be tuned with the environment variables `QC_MODE` (`auto` or
`window`), `THREADS`, `REF`, and `OUT`.