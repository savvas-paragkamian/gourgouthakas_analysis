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

cd ../Streptomyces_genomes_db

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

The files Streptomyces_gene_cluster_output.txt and BiobankStreptomyces_gene_cluster_output.txt were copied inside the "Streptomyces_genomes_anvio_Graphs". Then, the graphs were created in R (version 4.5.1) using the [pangenome_graphs.R](scripts/pangenome_graphs.R) script. This script was executed inside the "Streptomyces_genomes_anvio_Graphs" directory. 