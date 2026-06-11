# Gourgouthakas cave analysis

Here we analyse the data from the biobank of microbial isolates from 
the Gourgouthakas cave. 

## Genome assembly + QC

The cave *Streptomyces* genomes (`SRL*`) are assembled from short + long reads by
[`scripts/genomes/assembly.sh`](scripts/genomes/assembly.sh): `fastp`/`fastplong`
quality filtering, `unicycler` hybrid assembly run in parallel with GNU
`parallel`, and `quast` statistics. (The `bakta` annotation and `gtdb-tk`
taxonomy steps further down that script still run on the lab server.)

These tools are containerized in
[`scripts/genomes/Containerfile`](scripts/genomes/Containerfile) (one image, one
micromamba env per tool, since their pinned dependencies do not co-solve):

```
podman build -t gourgouthakas-genomes -f scripts/genomes/Containerfile .
# version smoke test:
podman run --rm gourgouthakas-genomes
# run the assembly + QC stage (mount the repo; tools are picked per env):
podman run --rm -v "$PWD":/work -w /work gourgouthakas-genomes \
    bash scripts/genomes/assembly.sh <dirs_file> <path>
```

Each tool is invoked as `micromamba run -n <env> <cmd>` (envs: `qc`,
`unicycler`, `autocycler`, `quast`); `assembly.sh` uses
`RUN=${RUN:-micromamba run -n}` so the prefix can be overridden off-container.

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

### Connecting SILVA and GTDB taxonomy

To relate the SILVA assignments to the [GTDB](https://gtdb.ecogenomic.org/)
taxonomy, the reads are additionally classified against the GTDB SSU
representative database (`bac120_ssu_reps` + `ar53_ssu_reps`, latest release)
using the same `vsearch --usearch_global` approach, and the two taxonomies are
joined per read. This is a standalone step run after the main pipeline:

```
podman run --rm -v "$PWD":/work -w /work sanger16s bash scripts/sanger/05_silva_gtdb.sh
```

It is implemented in [scripts/sanger/05_silva_gtdb.sh](scripts/sanger/05_silva_gtdb.sh)
(GTDB classification via [02_assign.py](scripts/sanger/02_assign.py), `--db gtdb`)
and [scripts/sanger/05_join.py](scripts/sanger/05_join.py) (join + concordance).
For every batch it writes to `results/sanger/<batch>/`:

```
<batch>.tax.gtdb.csv         GTDB taxonomy per read
<batch>.silva_vs_gtdb.csv    SILVA and GTDB lineage side by side per read
<batch>.concordance.csv      per-read genus/species agreement
```
plus a roll-up `results/sanger/silva_vs_gtdb.summary.csv`.

Because GTDB appends polyphyly suffixes (e.g. `Bacillus_A`, `Brevibacillus_B`),
the concordance compares genus names both strictly and "loosely" (suffix
stripped), and compares the species epithet ignoring SILVA placeholders such as
`sp.`/`uncultured`. Genuine GTDB reclassifications are surfaced this way — for
the cave isolates the largest is the split of the *Pseudomonas* complex into
*Aquipseudomonas*, *Stutzerimonas*, and *Neopusillimonas*.

### Combined taxonomy table

[scripts/sanger/08_taxonomy_table.py](scripts/sanger/08_taxonomy_table.py)
merges the SILVA, NCBI and GTDB assignments of every isolate into a single
`results/taxonomy_per_microbe.tsv`, with each database's reference identifier
(best hit), % identity, genus, species and full lineage side by side, and
**inner-joins** the cave-isolate metadata:

```
podman run --rm -v "$PWD":/work -w /work sanger16s \
  python scripts/sanger/08_taxonomy_table.py
```

It reads the **per-batch** isoTAX outputs rather than the merged ones, because
read ids restart per plate (e.g. `1_27f-A01` exists on several plates) and the
merged step keeps only the first of each collision. Each microbe is therefore
keyed by `(plate, microbe_id)`. Columns: `plate, microbe_id, stab,` then
`{silva,ncbi,gtdb}_{id,pct_id,genus,species,lineage}` (plus an `ncbi_taxid`
column, the NCBI accession mapped through `data/ref/ncbi_16S.acc2taxid.tsv`),
and finally every metadata column of `data/gourgouthakas-cave-isolates.csv`.

The metadata is joined by **`stab`** — the numeric isolate id, which is the
leading number of `microbe_id` (`1280_27f-F12` → `1280`, `1056-Premixed` →
`1056`). It is an **inner join**: only isolates present in the metadata sheet
are kept (the sheet covers stabs ≥ 475, so sequenced isolates absent from it
are dropped from the table and the depth abundances).

From the joined table the script also writes the per-depth abundance matrices
the tree figures consume, **one per taxonomy**:
`results/gourgouthakas_depth_table.gtdb.tsv` and `…silva.tsv`. Each has one row
per taxon (that taxonomy's genus and species labels) and one column per
Gourgouthakas sampling depth (`0, -39, -220, -418, -678, -713, -900, -1050,
-1100` m), counting the isolates of that taxon recovered at each depth. The CSV
`depth` column holds the positive magnitude; the table columns are the signed
depths. The two taxonomies are kept in **separate files** so a figure only sees
its own labels — otherwise the GTDB Pseudomonas-complex genera
(*Aquipseudomonas*, *Stutzerimonas*, …) would each also be counted under
*Pseudomonas* (their SILVA genus) in the GTDB figure, inflating and duplicating
them. **Run `08` before `06`/`07`.**

#### Whole-genome (gtdb-tk) taxonomy

For isolates that were also whole-genome sequenced and classified with
[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk), the genome-based assignment is
more reliable than the single 16S gene, so it takes precedence in the **GTDB**
outputs (the SILVA table is left untouched, since the WGS call is GTDB). The
GTDB-Tk ANI summaries are placed under `results/genomes/`
(`gtdbtk.ani_summary*.tsv`, one per assembly; the `user_genome` column is
`<stab>.unicycler_hybrid_assembly`). `08_taxonomy_table.py` takes the best
(highest `skani_ani`) reference hit of each assembly and, keyed by `stab`,
**overrides** that isolate's Sanger GTDB genus/species/lineage where the two
differ, and **inserts** any WGS-only isolate (one never recovered by Sanger 16S)
into both `taxonomy_per_microbe.tsv` (plate `WGS`) and the GTDB depth table. It
also writes `results/genomes/wgs_gtdb.tsv` (the per-isolate WGS GTDB call, with
the reference accession normalised to the GTDB master-tree prefix) for the tree
figure to consume.

### Truncated GTDB master tree

The GTDB **bac120 master tree** is pruned down to just the reference genomes the
isolates were assigned to (their GTDB best hits), so the isolates can be placed
in the context of the GTDB phylogeny. This follows the `keep.tip` / `ggtree`
approach of [scripts/taxonomy_tree.R](scripts/taxonomy_tree.R) and runs in R
inside the same container (which now also carries `r-tidyverse`, `r-ape`,
`ggtree`, `treeio`, `tidytree`, `ggtreeExtra`):

```
podman run --rm -v "$PWD":/work -w /work sanger16s \
  Rscript scripts/sanger/06_gtdb_tree.R [gtdb_tax.csv] [bac120.tree] [bac120_taxonomy.tsv] [out_prefix]
```

With no arguments it uses the merged GTDB results
(`results/sanger/merged/merged.tax.gtdb.csv`) and the master tree + taxonomy
downloaded into `data/ref/gtdb/` by [setup_db.sh](scripts/sanger/setup_db.sh).
[scripts/sanger/06_gtdb_tree.R](scripts/sanger/06_gtdb_tree.R) reads the GTDB
best-hit genome of each isolate, prunes the master tree with `keep.tip`, and
lays the tree out compressed into the left third of an **A4 page** next to a
`total` column (isolates per taxon) and a table of abundances across the
**9 Gourgouthakas sampling depths**
(0, -39, -220, -418, -678, -713, -900, -1050, -1100 m). Tip points are coloured
by phylum and sized by the number of isolates. **Figures are written to
`plots/`**; the pruned tree and tip table go next to the input:

```
plots/*.gtdb_master_tree[.by_genus].pdf / .png   A4 tree + depth-abundance table
plots/*.gtdb_master_tree[.by_genus].circular.png circular tree, coloured by genus
results/sanger/merged/*.nwk                       the pruned tree (Newick)
results/sanger/merged/*.tips.csv                  tip table (genome, lineage, n_isolates)
```

The depth abundances are read from the per-taxonomy tables
`results/gourgouthakas_depth_table.{gtdb,silva}.tsv` (taxon + 9 depth columns),
**generated by `08_taxonomy_table.py`** from the isolate metadata — run `08`
before `06`/`07`. Each figure reads the table matching the taxonomy it labels
with (`06` → `.gtdb`, `06` with `TAX=silva` and `07` → `.silva`). `06`/`07` only
read the table (they never write it): they keep only the taxa present in it
(dropping any with no depth metadata) and set each taxon's `total` column and
tip size to its depth-row sum, so the `total` equals the sum of the depth cells.
If the table is absent the scripts stop and tell you to run `08`. These files
are kept under version control even though `results/` is otherwise git-ignored.

When a whole-genome (GTDB-Tk) call is available (`results/genomes/wgs_gtdb.tsv`,
written by `08`), `06` **grafts** that isolate onto the pruned tree as a
synthetic hit on its GTDB best-hit genome (GTDB labelling only). Because the
bac120 master tree carries fewer genomes than GTDB as a whole, the exact ANI
reference is sometimes not a tip; `06` then falls back to a tree tip of the same
GTDB genus, so the isolate still appears under its genus (e.g. the cave
*Nocardiopsis* isolate, whose reference genome is absent from the tree).

Set `COLLAPSE=genus` for the compact one-tip-per-genus view (a representative
genome per GTDB genus, tip size = total isolates in the genus); outputs get a
`.by_genus` suffix so both views coexist. This is the version that fits one A4
page — for the merged set it collapses the 186 genomes to 69 genera:

```
podman run --rm -e COLLAPSE=genus -v "$PWD":/work -w /work sanger16s \
  Rscript scripts/sanger/06_gtdb_tree.R
```

Add `TAX=silva` to label/collapse this same GTDB-pruned tree by the **SILVA**
taxonomy instead (the tree stays the pruned GTDB master tree; genus/phylum come
from the SILVA per-read assignments, written as `*.by_genus.silva.*`):

```
podman run --rm -e COLLAPSE=genus -e TAX=silva -v "$PWD":/work -w /work sanger16s \
  Rscript scripts/sanger/06_gtdb_tree.R
```

The same genus tree + depth table can instead be drawn from the **de novo
FastTree** of the isolate cluster representatives (`04_tree.sh`, i.e. topology
from the isolates' own 16S sequences rather than the GTDB reference), using the
**SILVA** taxonomy for the genus labels, with
[scripts/sanger/07_fasttree_genus_tree.R](scripts/sanger/07_fasttree_genus_tree.R),
saved under the separate name `merged.fasttree_genus_tree.*` (pass
`merged.tax.gtdb.csv` as the first argument to use GTDB instead):

```
podman run --rm -v "$PWD":/work -w /work sanger16s \
  Rscript scripts/sanger/07_fasttree_genus_tree.R
```