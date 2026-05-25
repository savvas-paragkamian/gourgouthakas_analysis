######################################################################################################
# script name: pangenome_graphs.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to create the barplots and the UpSet plots of the pangenome analysis
######################################################################################################
# usage:./pangenome_graphs.R
# complete path: /home/nik_arapitsas/Documents/gourgouthakas_analysis/scripts/pangenome_graphs.R
######################################################################################################


# 1) Barplot for the Streptomyces pangenome

library(tidyverse)

# Input file
pangenome_file <- "Streptomyces_gene_cluster_output.txt"

# Read table
pangenome_table <- read_delim(pangenome_file, delim = "\t", col_types = cols())

# Count number of genomes in the pangenome analysis
n_genomes <- length(unique(pangenome_table$genome_name))

# Summarize gene clusters into core, accessory, and species-specific
pangenome_data <- pangenome_table %>%
  distinct(gene_cluster_id, genome_name) %>%
  count(gene_cluster_id) %>%
  mutate(category = case_when(
    n == n_genomes ~ "core",
    n == 1 ~ "species-specific",
    TRUE ~ "accessory"
  )) %>%
  count(category) %>%
  mutate(
    total = sum(n),
    percent = round(n / total * 100, 1),
    group = "Streptomyces pangenome"
  ) %>%
  select(group, category, n, percent)

# Set category order
pangenome_data$category <- factor(
  pangenome_data$category,
  levels = c("species-specific", "accessory", "core")
)

# Make one stacked bar
pangenome_barplot <- ggplot(
  pangenome_data,
  aes(x = "Streptomyces pangenome", y = percent, fill = category)
) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(
    aes(label = paste0(percent, "%")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 4.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    breaks = c("core", "accessory", "species-specific"),
    values = c(
      "core" = "#56B4E9",
      "accessory" = "#CC79A7",
      "species-specific" = "#E69F00"
    )
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100)
  ) +
  labs(
    x = "Streptomyces pangenome",
    y = "Pangenome composition (%)",
    fill = "Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),

    # remove the x-axis tick label
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),

    axis.text.y = element_text(
      size = 12,
      color = "black",
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 16,
      color = "black",
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(margin = margin(r = 12)),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save plot
ggsave(
  filename = "Streptomyces_pangenome_stacked_barplot.png",
  plot = pangenome_barplot,
  height = 20, 
  width = 25,
  dpi = 300, 
  units="cm",
  device="png"
)

# 2) UpSet Plot for Streptomyces pangenome

library(UpSetR)

# Create presence/absence matrix
Streptomyces_matrix <- pangenome_table %>%
  distinct(gene_cluster_id, genome_name) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = genome_name, values_from = present, values_fill = 0) %>%
  column_to_rownames("gene_cluster_id")

colnames(Streptomyces_matrix) <- gsub("_(\\d+)$", ".\\1", colnames(Streptomyces_matrix))

# Map accessions to readable isolate names
name_map <- c(
  "GCF_000009765.2" = "S. avermitilis MA-4680",
  "GCF_000203835.1" = "S. coelicolor A3(2)",
  "GCF_003248355.1" = "Streptomyces sp. AC1-42W",
  "GCF_003287935.1" = "Streptomyces sp. ICC1",
  "GCF_003287915.1" = "Streptomyces sp. ICC4",
  "GCF_003054555.1" = "S. lunaelactis MM109",
  "GCF_048572725.1" = "S. boninensis TBRC 7755",
  "GCF_049310085.1" = "S. tritrimontium 4.24",
  "GCF_049310105.1" = "S. tritrimontium 2.9",
  "GCF_049310065.1" = "S. tatrensis 5.8",
  "GCF_003248335.1" = "Streptomyces sp. AC1-42T",
  "SRL740" = "★ Streptomyces sp. SRL740",
  "SRL742" = "★ Streptomyces sp. SRL742",
  "SRL1060" = "★ Streptomyces sp. SRL1060"
)

colnames(Streptomyces_matrix) <- ifelse(
  colnames(Streptomyces_matrix) %in% names(name_map),
  name_map[colnames(Streptomyces_matrix)],
  colnames(Streptomyces_matrix)
)

# UpSet plot
png("Streptomyces_upset.png", width = 4900, height = 4000, res = 300)

upset(
  Streptomyces_matrix,
  sets = colnames(Streptomyces_matrix),
  order.by = "freq",
  mainbar.y.label = "\n\n\n\n\n\n\n\nGene Cluster Intersections - Streptomyces",
  sets.x.label = "Gene Clusters per Genome",
  text.scale = c(1.3, 1.5, 1.5, 1.5, 1.5, 1.3)
)

dev.off()


# 3) Barplot for the biobank's Streptomyces pangenome

biobank_pangenome_file <- "BiobankStreptomyces_gene_cluster_output.txt"

# Read table
biobank_pangenome_table <- read_delim(biobank_pangenome_file, delim = "\t", col_types = cols())

# Count number of genomes in the pangenome analysis
n_genomes <- length(unique(biobank_pangenome_table$genome_name))

# Summarize gene clusters into core, accessory, and species-specific
biobank_pangenome_data <- biobank_pangenome_table %>%
  distinct(gene_cluster_id, genome_name) %>%
  count(gene_cluster_id) %>%
  mutate(category = case_when(
    n == n_genomes ~ "core",
    n == 1 ~ "species-specific",
    TRUE ~ "accessory"
  )) %>%
  count(category) %>%
  mutate(
    total = sum(n),
    percent = round(n / total * 100, 1),
    group = "Biobank Streptomyces pangenome"
  ) %>%
  select(group, category, n, percent)

# Set category order
biobank_pangenome_data$category <- factor(
  biobank_pangenome_data$category,
  levels = c("species-specific", "accessory", "core")
)

# Make one stacked bar
biobank_pangenome_barplot <- ggplot(
  biobank_pangenome_data,
  aes(x = "Biobank Streptomyces pangenome", y = percent, fill = category)
) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(
    aes(label = paste0(percent, "%")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 4.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    breaks = c("core", "accessory", "species-specific"),
    values = c(
      "core" = "#56B4E9",
      "accessory" = "#CC79A7",
      "species-specific" = "#E69F00"
    )
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 100)
  ) +
  labs(
    x = "Biobank Streptomyces pangenome",
    y = "Pangenome composition (%)",
    fill = "Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),

    # remove the duplicated x-axis tick label
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),

    axis.text.y = element_text(
      size = 12,
      color = "black",
      family = "sans"
    ),
    axis.title.x = element_text(
      size = 16,
      color = "black",
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(margin = margin(r = 12)),
    plot.margin = margin(10, 10, 10, 10)
  )
# Save plot
ggsave(
  filename = "Biobank_Streptomyces_pangenome_stacked_barplot.png",
  plot = biobank_pangenome_barplot,
  height = 20, 
  width = 25,
  dpi = 300, 
  units="cm",
  device="png"
)

# 4) UpSet Plot for the biobank's Streptomyces pangenome

# Create presence/absence matrix
biobank_Streptomyces_matrix <- biobank_pangenome_table %>%
  distinct(gene_cluster_id, genome_name) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = genome_name, values_from = present, values_fill = 0) %>%
  column_to_rownames("gene_cluster_id")

biobank_name_map <- c(
  "SRL740" = "Streptomyces sp. SRL740",
  "SRL742" = "Streptomyces sp. SRL742",
  "SRL1060" = "Streptomyces sp. SRL1060"
)

colnames(biobank_Streptomyces_matrix) <- ifelse(
  colnames(biobank_Streptomyces_matrix) %in% names(biobank_name_map),
  biobank_name_map[colnames(biobank_Streptomyces_matrix)],
  colnames(biobank_Streptomyces_matrix)
)

# UpSet plot
png("Biobank_Streptomyces_upset.png", width = 4900, height = 4000, res = 300)

upset(
  biobank_Streptomyces_matrix,
  sets = colnames(biobank_Streptomyces_matrix),
  order.by = "freq",
  mainbar.y.label = "\n\n\n\n\n\nGene Cluster Intersections - Biobank Streptomyces",
  sets.x.label = "Gene Clusters per Genome",
  text.scale = c(1.2, 1.2, 1, 0.9, 1.0, 1.1)
)

dev.off()