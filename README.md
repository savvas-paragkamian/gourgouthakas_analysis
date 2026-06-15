# Companion analysis of Gourgouthakas cave biobank

Biobank of microbial isolates from the **Gourgouthakas cave** (Crete, Greece),
contains more than 820 isolates from samples across nine depths down to −1100 m.
This repository holds the
reproducible analysis of that biobank, in three tracks:

1. **Genome assembly, QC and Annotation** of the cave *Streptomyces* , *Nocardiopsis* and *Pseudomonas* isolates (hybrid
   short + long read assembly).
2. **Pangenomics** (anvi'o) of those *Streptomyces*, against selected public
   genomes and among themselves.
3. **Taxonomic identification** of the Sanger 16S–sequenced isolates (an
   [isolateR](https://github.com/bdaisley/isolateR) style pipeline), their placement on the GTDB phylogeny, and
   per-depth abundance figures.

Everything runs from **podman** containers with the logic in `scripts/`, so the
analyses reproduce on machines where the bioinformatics tools won't install
natively. **Step-by-step commands and reference detail live in
[`notes.md`](notes.md);** this README is the overview.

## Pipelines

| Track | Scripts | Container |
|---|---|---|
| Genome assembly + QC | [`scripts/genomes/assembly.sh`](scripts/genomes/assembly.sh) | `gourgouthakas-genomes` |
| Pangenomics (anvi'o) | [`pangenome_graphs.R`](scripts/genomes/pangenome_graphs.R), [`pangenome_functions.R`](scripts/genomes/pangenome_functions.R) | anvi'o (lab server) |
| Sanger 16S taxonomy + GTDB figures | [`scripts/sanger/`](scripts/sanger) (`01`–`08`) | `sanger16s` |

## Genome assembly + QC

The cave *Streptomyces* genomes (`SRL*`) are assembled from short + long reads by
[`scripts/genomes/assembly.sh`](scripts/genomes/assembly.sh): `fastp`/`fastplong`
quality filtering, `unicycler` hybrid assembly (GNU `parallel`), and `quast`
statistics. (The `bakta` annotation and `gtdb-tk` taxonomy steps further down
that script still run on the lab server.) Tools are containerized in
[`scripts/genomes/Containerfile`](scripts/genomes/Containerfile) — one image,
one micromamba env per tool, since their pinned dependencies don't co-solve.

```
podman build -t gourgouthakas-genomes -f scripts/genomes/Containerfile .
podman run --rm -v "$PWD":/work -w /work gourgouthakas-genomes \
    bash scripts/genomes/assembly.sh <dirs_file> <path>
```

→ [Full commands in `notes.md`](notes.md#genome-assembly--qc).

## Pangenome analysis (anvi'o)

Two anvi'o pangenomes are built: the three cave *Streptomyces* isolates against
eleven selected public *Streptomyces* genomes, and the three isolates alone.
Each genome gets an anvi'o contigs database annotated with HMMs, NCBI COGs, KEGG
KOfams, Pfams and SCG taxonomy; `anvi-pan-genome` (`--minbit 0.5
--mcl-inflation 10`) computes the gene clusters. Gene-cluster presence/absence
and per-isolate functions are exported and analysed in R with
[`pangenome_graphs.R`](scripts/genomes/pangenome_graphs.R) and
[`pangenome_functions.R`](scripts/genomes/pangenome_functions.R).

→ [Full anvi'o command walkthrough in `notes.md`](notes.md#pangenome-analysis-with-anvio).

## Taxonomic identification from Sanger 16S

Isolates were identified by Sanger 16S (27f primer); raw `.ab1` traces live under
`data/` (one directory per plate, plus the premixed set). The
[isolateR](https://github.com/bdaisley/isolateR) workflow (Daisley et al. 2024,
*Bioinformatics* 40(7):btae448) was rebuilt with vibe coding and manual evaluation
as command-line steps in the
`sanger16s` podman container, because the R package would not install:

| isolateR step | Script | Does |
|---|---|---|
| `isoQC` | [`01_isoqc.py`](scripts/sanger/01_isoqc.py) | quality-trim `.ab1` → FASTA + PASS/FAIL |
| `isoTAX` | [`02_isotax.sh`](scripts/sanger/02_isotax.sh) | global-alignment taxonomy vs NCBI 16S + SILVA |
| `isoLIB` | [`03_isolib.sh`](scripts/sanger/03_isolib.sh) | 99.5 % strain library |
| tree | [`04_tree.sh`](scripts/sanger/04_tree.sh) | mafft + FastTree of representatives |

```
podman build -t sanger16s -f scripts/sanger/Containerfile .
podman run --rm -v "$PWD":/work -w /work sanger16s bash scripts/sanger/run.sh
```

[`run.sh`](scripts/sanger/run.sh) processes each plate as a batch and a merged
cross-plate set, writing to `results/sanger/<batch>/`. Reference DBs (NCBI 16S,
SILVA SSU NR99, GTDB) are fetched once by
[`setup_db.sh`](scripts/sanger/setup_db.sh).

→ [Outputs, parameters, and the isolateR-mirrored cutoffs in `notes.md`](notes.md#sanger-16s-pipeline).

### GTDB placement and per-depth figures

Downstream of the per-plate taxonomy:

- **`05`** ([`05_silva_gtdb.sh`](scripts/sanger/05_silva_gtdb.sh)) also classifies
  each read against GTDB and reports per-read SILVA↔GTDB concordance.
- **`08`** ([`08_taxonomy_table.py`](scripts/sanger/08_taxonomy_table.py)) merges
  the SILVA/NCBI/GTDB calls per isolate into `results/taxonomy_per_microbe.tsv`
  (inner-joined to the cave metadata) and writes the per-depth abundance tables
  the figures consume (`results/gourgouthakas_depth_table.{gtdb,silva}.tsv`, one
  per taxonomy). Where an isolate was **whole-genome sequenced** and classified
  with GTDB-Tk (`results/genomes/gtdbtk.ani_summary*.tsv`), that genome-based
  call **overrides** the 16S GTDB call (SILVA is left untouched).
- **`06`/`07`** ([`06_gtdb_tree.R`](scripts/sanger/06_gtdb_tree.R),
  [`07_fasttree_genus_tree.R`](scripts/sanger/07_fasttree_genus_tree.R)) draw the
  final figure: the GTDB bac120 master tree (or the de novo FastTree) pruned to
  the isolates' genera, beside a heatmap of abundance across the nine cave
  sampling depths. WGS isolates are grafted onto the tree. Figures → `plots/`.

```
podman run --rm -v "$PWD":/work -w /work sanger16s python scripts/sanger/08_taxonomy_table.py
podman run --rm -e COLLAPSE=genus -v "$PWD":/work -w /work sanger16s Rscript scripts/sanger/06_gtdb_tree.R
```

→ [Concordance internals, table columns, the WGS override, and the figure
variants (`TAX=silva`, `COLLAPSE`, `07`) in `notes.md`](notes.md#gtdb-placement-and-depth-figures).

## Phytopathogen antagonism

[`scripts/pathogen_inhibition.R`](scripts/pathogen_inhibition.R) (run in the
`sanger16s` image) has two parts:

- **Part 1 — in vitro screen** of the isolates against six phytopathogens
  (`data/in_vitro_phytopathogens_inhibition.txt`): inhibition-score heatmap and
  per-genus bar/bubble figures → `plots/`.
- **Part 2 — ex-vivo *Botrytis cinerea* biocontrol** (`data/ex-vivo-inhibition_
  B.c._SRL917_gourgouthakas.xlsx`, sheet `raw`; four treatments × 20 reps). Per
  d.p.i. **one-way ANOVA + Tukey HSD** (each treatment vs *B. cinerea* alone)
  and an `audpc2()` trapezoidal **AUDPC** (anchored at inoculation day 0,
  reproducing the sheet exactly). Stats → `results/ex_vivo_anova.tsv`,
  `results/ex_vivo_tukey.tsv`; grouped per-d.p.i. and AUDPC bar plots (SE error
  bars, Tukey significance stars) → `plots/ex_vivo_barplot_{dpi,audpc}.png`.

**SRL917** suppresses disease at every time point (`***`, AUDPC 24 vs 77 for the
pathogen alone); the **X** product is significant only early/mid and by AUDPC,
but not at 5–6 d.p.i. (small effect swamped by the pathogen group's high
late-stage variance).

```
podman run --rm -v "$PWD":/work -w /work/scripts sanger16s Rscript pathogen_inhibition.R
```
