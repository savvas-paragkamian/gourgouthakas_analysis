library(tidyverse)

# ============================================================
# 1. Read exported gene cluster table
# ============================================================

gc_table <- read_tsv("Streptomyces_gene_clusters.tsv", col_types = cols())

# ============================================================
# 2. Define the biobank's cave isolates
# ============================================================

isolates <- c("SRL1060", "SRL740", "SRL742")

# ============================================================
# 3. Presence/absence per gene cluster
# ============================================================

cluster_presence <- gc_table %>%
  distinct(gene_cluster_name, genome_name) %>%
  group_by(gene_cluster_name) %>%
  summarise(
    genomes_present = paste(sort(unique(genome_name)), collapse = ";"),
    n_genomes = n_distinct(genome_name),
    present_in_SRL1060 = "SRL1060" %in% genome_name,
    present_in_SRL740  = "SRL740" %in% genome_name,
    present_in_SRL742  = "SRL742" %in% genome_name,
    .groups = "drop"
  )

# ============================================================
# 4. Gene clusters unique to SRL1060 only
# ============================================================

SRL1060_unique_clusters <- cluster_presence %>%
  filter(
    n_genomes == 1,
    present_in_SRL1060
  )

# ============================================================
# 5. Gene clusters present only in SRL740, SRL742, and SRL1060
# ============================================================

SRL740_742_1060_shared_unique_clusters <- cluster_presence %>%
  filter(
    n_genomes == 3,
    present_in_SRL1060,
    present_in_SRL740,
    present_in_SRL742
  )

# ============================================================
# 6. Gene clusters present only in SRL740 and SRL742
# ============================================================

SRL740_742_shared_unique_clusters <- cluster_presence %>%
  filter(
    n_genomes == 2,
    present_in_SRL740,
    present_in_SRL742
  )

# ============================================================
# 7. Extract genes belonging to these clusters
# ============================================================

SRL1060_unique_genes <- gc_table %>%
  filter(gene_cluster_name %in% SRL1060_unique_clusters$gene_cluster_name)

SRL740_742_1060_shared_unique_genes <- gc_table %>%
  filter(gene_cluster_name %in% SRL740_742_1060_shared_unique_clusters$gene_cluster_name)

SRL740_742_shared_unique_genes <- gc_table %>%
  filter(gene_cluster_name %in% SRL740_742_shared_unique_clusters$gene_cluster_name)

# ============================================================
# 8. Save gene cluster and gene lists
# ============================================================

write_tsv(
  SRL1060_unique_clusters,
  "SRL1060_unique_gene_clusters.tsv"
)

write_tsv(
  SRL740_742_1060_shared_unique_clusters,
  "SRL740_SRL742_SRL1060_shared_unique_gene_clusters.tsv"
)

write_tsv(
  SRL740_742_shared_unique_clusters,
  "SRL740_SRL742_shared_unique_gene_clusters.tsv"
)

write_tsv(
  SRL1060_unique_genes,
  "SRL1060_unique_genes.tsv"
)

write_tsv(
  SRL740_742_1060_shared_unique_genes,
  "SRL740_SRL742_SRL1060_shared_unique_genes.tsv"
)

write_tsv(
  SRL740_742_shared_unique_genes,
  "SRL740_SRL742_shared_unique_genes.tsv"
)

# ============================================================
# 9. Read functional annotation files for the three isolates
# ============================================================

SRL1060_functions <- read_tsv("SRL1060_functions.tsv", col_types = cols()) %>%
  mutate(genome_name = "SRL1060")

SRL740_functions <- read_tsv("SRL740_functions.tsv", col_types = cols()) %>%
  mutate(genome_name = "SRL740")

SRL742_functions <- read_tsv("SRL742_functions.tsv", col_types = cols()) %>%
  mutate(genome_name = "SRL742")


# Combine all functions
functions <- bind_rows(
  SRL1060_functions,
  SRL740_functions,
  SRL742_functions
)


# Rename gene_callers_id to match gene_caller_id in gc table
functions <- functions %>%
  rename(gene_caller_id = gene_callers_id)


# ============================================================
# 10. Join functional annotations
# ============================================================

SRL1060_unique_genes_annotated <- SRL1060_unique_genes %>%
  left_join(
    functions,
    by = c("genome_name", "gene_caller_id")
  )

SRL740_742_1060_shared_unique_genes_annotated <- SRL740_742_1060_shared_unique_genes %>%
  left_join(
    functions,
    by = c("genome_name", "gene_caller_id")
  )

SRL740_742_shared_unique_genes_annotated <- SRL740_742_shared_unique_genes %>%
  left_join(
    functions,
    by = c("genome_name", "gene_caller_id")
  )


# ============================================================
# 11. Save annotated gene files
# ============================================================

write_tsv(
  SRL1060_unique_genes_annotated,
  "SRL1060_unique_genes_annotated.tsv"
)

write_tsv(
  SRL740_742_1060_shared_unique_genes_annotated,
  "SRL740_SRL742_SRL1060_shared_unique_genes_annotated.tsv"
)

write_tsv(
  SRL740_742_shared_unique_genes_annotated,
  "SRL740_SRL742_shared_unique_genes_annotated.tsv"
)

# ============================================================
# 12. Functional summaries
# ============================================================

SRL1060_unique_genes_annotated %>%
  filter(!is.na(`function`)) %>%
  count(source, `function`, sort = TRUE) %>%
  write_tsv("SRL1060_unique_function_summary.tsv")

SRL740_742_1060_shared_unique_genes_annotated %>%
  filter(!is.na(`function`)) %>%
  count(source, `function`, sort = TRUE) %>%
  write_tsv("SRL740_SRL742_SRL1060_shared_unique_function_summary.tsv")

SRL740_742_shared_unique_genes_annotated %>%
  filter(!is.na(`function`)) %>%
  count(source, `function`, sort = TRUE) %>%
  write_tsv("SRL740_SRL742_shared_unique_function_summary.tsv")

# View top 30 functions
SRL1060_unique_genes_annotated %>%
  filter(!is.na(`function`)) %>%
  count(`function`, sort = TRUE) %>%
  head(30)

SRL740_742_1060_shared_unique_genes_annotated %>%
  filter(!is.na(`function`)) %>%
  count(`function`, sort = TRUE) %>%
  head(30)

SRL740_742_shared_unique_genes_annotated %>%
  filter(!is.na(`function`)) %>%
  count(`function`, sort = TRUE) %>%
  head(30)
