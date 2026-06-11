#!/usr/bin/env Rscript

###############################################################################
# script name: 06_gtdb_tree.R
# developed by: Savvas Paragkamian (Sanger 16S pipeline)
# framework: SarrisLab
###############################################################################
# GOAL:
# Truncate the GTDB bac120 master tree to the reference genomes that the cave
# isolates were assigned to (GTDB best hits from isoTAX-vs-GTDB) and plot it.
# Follows the keep.tip / ggtree approach of scripts/taxonomy_tree.R.
###############################################################################
# usage (inside the sanger16s container, from repo root):
#   Rscript scripts/sanger/06_gtdb_tree.R \
#       [gtdb_tax.csv] [bac120.tree] [bac120_taxonomy.tsv] [out_prefix] [depth.tsv]
#   COLLAPSE=genus  -> one representative tip per genus (compact, one A4 page)
#   TAX=silva       -> label/collapse by SILVA taxonomy (tree stays the pruned
#                      GTDB master tree); default labels are GTDB.
# Figures go to plots/; the depth table is created (random placeholder) in
# results/ if missing -- edit it with real abundances and re-run.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggtree)
  library(tidytree)
  library(treeio)
})

args <- commandArgs(trailingOnly = TRUE)
tax_csv   <- ifelse(length(args) >= 1, args[1],
                    "results/sanger/merged/merged.tax.gtdb.csv")
tree_file <- ifelse(length(args) >= 2, args[2], "data/ref/gtdb/bac120.tree")
gtdb_tax  <- ifelse(length(args) >= 3, args[3],
                    "data/ref/gtdb/bac120_taxonomy.tsv")
out_pref  <- ifelse(length(args) >= 4, args[4],
                    file.path(dirname(tax_csv),
                              sub("\\.tax\\.gtdb\\.csv$", "", basename(tax_csv)) %>%
                                paste0(".gtdb_master_tree")))

# COLLAPSE=genus keeps one representative tip per GTDB genus (compact tree);
# default "none" keeps one tip per best-hit genome.
collapse <- tolower(Sys.getenv("COLLAPSE", "none"))
if (collapse == "genus") out_pref <- paste0(out_pref, ".by_genus")

# TAX=silva labels the GTDB-pruned tree with SILVA genus/phylum (the tree is
# still the GTDB master tree pruned to the GTDB best-hit genomes); default uses
# the GTDB reference taxonomy of those genomes.
tax_mode  <- tolower(Sys.getenv("TAX", "gtdb"))
silva_csv <- sub("\\.tax\\.gtdb\\.csv$", ".tax.silva.csv", tax_csv)
if (tax_mode == "silva") out_pref <- paste0(out_pref, ".silva")

message("input    : ", tax_csv)
message("tree     : ", tree_file)
message("taxonomy : ", tax_mode)
message("collapse : ", collapse)
message("out      : ", out_pref, ".{pdf,png}")

# ---- 1. our isolates' GTDB best-hit genomes --------------------------------
hits <- read_csv(tax_csv, show_col_types = FALSE) %>%
  filter(!is.na(best_hit), best_hit != "")

# one row per reference genome = a tip to keep, with how many isolates hit it
mode_chr <- function(x) {                 # most frequent non-empty value, or NA
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else names(which.max(table(x)))
}

per_genome <- hits %>%
  group_by(genome_id = best_hit) %>%
  summarise(n_isolates = n(),
            our_genus  = mode_chr(genus),
            mean_pct   = round(mean(suppressWarnings(as.numeric(pct_id)),
                                     na.rm = TRUE), 1),
            .groups = "drop")
message("distinct best-hit genomes: ", nrow(per_genome))

# ---- 2. keep only genomes present in the master tree, attach GTDB lineage ---
tree <- ape::read.tree(tree_file)
per_genome <- per_genome %>% filter(genome_id %in% tree$tip.label)
message("genomes present in master tree: ", nrow(per_genome))

ranks <- read_tsv(gtdb_tax, col_names = c("genome_id", "taxonomy"),
                  show_col_types = FALSE) %>%
  filter(genome_id %in% per_genome$genome_id) %>%
  mutate(phylum  = str_remove(str_extract(taxonomy, "p__[^;]+"), "p__"),
         class   = str_remove(str_extract(taxonomy, "c__[^;]+"), "c__"),
         genus   = str_remove(str_extract(taxonomy, "g__[^;]+"), "g__"),
         species = str_remove(str_extract(taxonomy, "s__[^;]+"), "s__")) %>%
  mutate(species = if_else(is.na(species) | species == "", genome_id, species))

genomes <- per_genome %>% left_join(ranks, by = "genome_id")

if (tax_mode == "silva") {
  # replace the GTDB genus/phylum of each genome with the SILVA assignment of
  # the reads that hit it (most frequent value across those reads)
  silva <- read_csv(silva_csv, show_col_types = FALSE)
  g2s <- hits %>% transmute(query, genome_id = best_hit) %>%
    inner_join(silva %>% transmute(query, sgenus = na_if(genus, ""),
                                   sphylum = na_if(phylum, "")), by = "query") %>%
    group_by(genome_id) %>%
    summarise(genus = mode_chr(sgenus), phylum = mode_chr(sphylum),
              .groups = "drop")
  genomes <- genomes %>% select(-genus, -phylum) %>%
    left_join(g2s, by = "genome_id")
  message("relabelled with SILVA taxonomy: ",
          sum(!is.na(genomes$genus)), "/", nrow(genomes), " genomes")
}

# ---- 3. choose which tips to keep (full vs one-per-genus) -------------------
if (collapse == "genus") {
  named   <- genomes %>% filter(!is.na(genus) & genus != "")
  unnamed <- genomes %>% filter(is.na(genus) | genus == "")   # keep individually

  # per-genus totals, and the representative genome (most isolates) per genus
  genus_tot <- named %>% group_by(genus) %>%
    summarise(n_isolates = sum(n_isolates), n_genomes = n_distinct(genome_id),
              phylum = first(phylum), class = first(class), .groups = "drop")
  reps <- named %>% arrange(desc(n_isolates)) %>% group_by(genus) %>%
    slice(1) %>% ungroup() %>% select(genome_id, genus) %>%
    left_join(genus_tot, by = "genus") %>%
    mutate(species = genus,
           tip_lab   = paste0(genus, "  (n=", n_isolates, ", ",
                              n_genomes, " genomes)"),
           short_lab = genus)

  un <- unnamed %>%
    mutate(n_genomes = 1L,
           tip_lab   = paste0(species, "  (n=", n_isolates, ")"),
           short_lab = species) %>%
    select(genome_id, genus, phylum, class, species,
           n_isolates, n_genomes, tip_lab, short_lab)

  chosen <- bind_rows(reps %>% select(all_of(names(un))), un)
} else {
  chosen <- genomes %>%
    mutate(n_genomes = 1L,
           tip_lab   = paste0(species, "  (n=", n_isolates, ")"),
           short_lab = species) %>%
    select(genome_id, genus, phylum, class, species,
           n_isolates, n_genomes, tip_lab, short_lab)
}
# ---- depth-abundance table (built by 08): keep only taxa it contains --------
# Real per-taxon per-depth counts live in results/gourgouthakas_depth_table.tsv
# (generated by 08 from the isolate metadata; 0 where a taxon is absent at a
# depth). Keep only figure taxa that have at least one isolate in it; drop the
# rest (no depth metadata) rather than fabricating counts.
depths    <- c(0, -39, -220, -418, -678, -713, -900, -1050, -1100)
# per-taxonomy depth table (08 writes .gtdb.tsv and .silva.tsv); label with the
# taxonomy this figure uses so genus counts match the tip labels.
depth_tsv <- ifelse(length(args) >= 5, args[5],
                    paste0("results/gourgouthakas_depth_table.", tax_mode, ".tsv"))
if (!file.exists(depth_tsv))
  stop("depth table not found: ", depth_tsv, " -- run 08_taxonomy_table.py first")
depth_df <- read_tsv(depth_tsv, show_col_types = FALSE)

dropped <- setdiff(chosen$short_lab, depth_df$taxon)
chosen  <- chosen %>% filter(short_lab %in% depth_df$taxon)
if (length(dropped) > 0)
  message("dropped ", length(dropped), " taxa absent from ",
          basename(depth_tsv), " (no depth metadata): ",
          paste(head(dropped, 6), collapse = ", "),
          if (length(dropped) > 6) ", ..." else "")

# Count only the supplied-table isolates: the "total" column and the tip size
# are the depth-table row sum, so total == sum of the depth cells (inner join).
depth_tot <- depth_df %>%
  transmute(short_lab = taxon,
            n_isolates = rowSums(across(all_of(as.character(depths)))))
chosen <- chosen %>% select(-n_isolates) %>%
  left_join(depth_tot, by = "short_lab")

message("tips kept: ", nrow(chosen),
        if (collapse == "genus") " (one per genus)" else " (one per genome)")

keep <- chosen$genome_id
if (length(keep) < 3) stop("Fewer than 3 matching tips; nothing to plot.")
pruned <- keep.tip(tree, keep)

tip_data <- tibble(label = pruned$tip.label) %>%
  left_join(chosen, by = c("label" = "genome_id"))

# ---- 5. tree (left, compressed) + table (right) on one A4 page -------------
fs_name <- 2.9; fs_cell <- 2.7; fs_hdr <- 3.4   # enlarged fonts for A4

# Okabe-Ito colour-blind-safe palette for the discrete phylum colours, ordered
# for maximum divergence between the few phyla shown (no two blues adjacent) and
# excluding the blue (#0072B2) used by the abundance heatmap below.
cb_palette <- c("#D55E00", "#009E73", "#E69F00", "#CC79A7",
                "#56B4E9", "#000000", "#999999")

p <- ggtree(pruned) %<+% tip_data +
  geom_tippoint(aes(color = phylum, size = n_isolates))

dd <- p$data %>% filter(isTip)
xr <- max(dd$x); n <- nrow(dd)

# layout in tree x-units: tree in the left ~fifth, then names, a "total" column,
# then the 9 depth columns. Wide gaps push the tree left (shrink the clades).
x_name <- xr * 1.05
dx     <- xr * 0.32          # column width
x_tot  <- xr * 2.05          # "total isolates" column
x0     <- x_tot + dx         # centre of the first depth column
hdr_y  <- n + 1

tot <- dd %>% distinct(short_lab, y, n_isolates) %>% mutate(xcell = x_tot)

long <- depth_df %>%
  pivot_longer(-taxon, names_to = "depth", values_to = "value") %>%
  mutate(depth = factor(depth, levels = as.character(depths)),
         xcell = x0 + (as.integer(depth) - 1) * dx) %>%
  inner_join(dd %>% select(short_lab, y), by = c("taxon" = "short_lab"))
# white->blue gradient is light at low values, dark at high -> white text only
# on the dark (high-abundance) cells, black on the light ones.
long <- long %>%
  mutate(txtcol = if_else(value > 0.55 * max(value, na.rm = TRUE),
                          "white", "black"))
hdr <- tibble(depth = as.character(depths),
              xcell = x0 + (seq_along(depths) - 1) * dx)
x_right <- max(hdr$xcell) + dx

p_tab <- p +
  geom_text(data = dd, aes(x = x_name, y = y, label = short_lab, color = phylum),
            hjust = 0, size = fs_name, fontface = "italic", family = "Ubuntu",
            inherit.aes = FALSE) +
  # total isolates per taxon (neutral fill to set it apart from the heatmap)
  geom_tile(data = tot, aes(x = xcell, y = y), fill = "grey92",
            width = dx * 0.92, height = 0.9, color = "grey80",
            linewidth = 0.1, inherit.aes = FALSE) +
  geom_text(data = tot, aes(x = xcell, y = y, label = n_isolates),
            size = fs_cell, fontface = "bold", inherit.aes = FALSE) +
  # depth-abundance heatmap columns
  geom_tile(data = long, aes(x = xcell, y = y, fill = value),
            width = dx * 0.92, height = 0.9, color = "grey80",
            linewidth = 0.1, inherit.aes = FALSE) +
  geom_text(data = long, aes(x = xcell, y = y, label = value),
            color = long$txtcol, size = fs_cell, inherit.aes = FALSE) +
  # column headers
  geom_text(data = hdr, aes(x = xcell, y = hdr_y, label = depth),
            angle = 90, hjust = 0, size = fs_hdr, fontface = "bold",
            inherit.aes = FALSE) +
  annotate("text", x = x_name, y = hdr_y, label = "taxon", hjust = 0,
           fontface = "bold", size = fs_hdr) +
  annotate("text", x = x_tot, y = hdr_y, label = "total", angle = 90,
           hjust = 0, fontface = "bold", size = fs_hdr) +
  scale_fill_gradient(low = "white", high = "#0072B2", name = "abundance") +
  scale_color_manual(values = cb_palette, name = "phylum",
                     na.value = "grey50") +
  scale_size_continuous(range = c(1, 5), name = "isolates") +
  coord_cartesian(xlim = c(0, x_right), ylim = c(0.5, n + 6), clip = "off") +
  theme_tree() +
  guides(
    fill  = guide_colorbar(title.position = "top", barheight = unit(0.3, "cm"),
                           barwidth = unit(3, "cm"), order = 1),
    color = guide_legend(title.position = "top", nrow = 2,
                         override.aes = list(size = 2.5), order = 2),
    size  = guide_legend(title.position = "top", nrow = 1, order = 3)) +
  # compact horizontal legend along the bottom so it fits within the page
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "top",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.margin = margin(1, 1, 1, 1),
        plot.margin = margin(2, 6, 2, 4))

# ---- 6. write figures to plots/, data next to the input --------------------
dir.create("plots", showWarnings = FALSE)
fig_pref <- file.path("plots", basename(out_pref))

ggsave(paste0(fig_pref, ".pdf"), p_tab, device = cairo_pdf,
       width = 21, height = 29.7, units = "cm", limitsize = FALSE)
ggsave(paste0(fig_pref, ".png"), p_tab, device = "png", dpi = 300,
       width = 21, height = 29.7, units = "cm", limitsize = FALSE)
# publication-quality 600 dpi TIFF (LZW-compressed) of the A4 tree + table
ggsave(paste0(fig_pref, ".tiff"), p_tab, device = "tiff", type = "cairo",
       dpi = 600, compression = "lzw",
       width = 21, height = 29.7, units = "cm", limitsize = FALSE)

# circular tree (no table), also to plots/
p_circ <- ggtree(pruned, layout = "circular") %<+% tip_data +
  geom_tippoint(aes(color = genus, size = n_isolates)) +
  geom_tiplab(aes(label = short_lab), size = 2.2, align = TRUE,
              linesize = 0.1, offset = 0.02, fontface = "italic",
              family = "Ubuntu") +
  scale_color_viridis_d(name = "genus") +
  scale_size_continuous(range = c(1, 6)) +
  theme(legend.position = "right",
        plot.margin = margin(2, 2, 2, 2, "cm"))
ggsave(paste0(fig_pref, ".circular.png"), p_circ, device = "png", dpi = 300,
       height = 40, width = 45, units = "cm", limitsize = FALSE)

write.tree(pruned, file = paste0(out_pref, ".nwk"))
write_csv(tip_data, paste0(out_pref, ".tips.csv"))

message("done: figures -> ", fig_pref, ".{pdf,png,tiff,circular.png}")
message("      data    -> ", out_pref, ".{nwk,tips.csv};  table -> ", depth_tsv)

