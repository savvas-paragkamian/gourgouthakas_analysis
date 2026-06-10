#!/usr/bin/env Rscript

###############################################################################
# script name: 07_fasttree_genus_tree.R
# developed by: Savvas Paragkamian (Sanger 16S pipeline)
# framework: SarrisLab
###############################################################################
# GOAL:
# Same genus-level tree + sampling-depth table as 06_gtdb_tree.R, but built from
# the *de novo* FastTree of the isolate cluster representatives (04_tree.sh)
# instead of the pruned GTDB master tree. Tips are collapsed to one
# representative per genus using the SILVA taxonomy. Saved under a separate name.
###############################################################################
# usage (inside the sanger16s container, from repo root):
#   Rscript scripts/sanger/07_fasttree_genus_tree.R \
#       [tax.csv] [reps.nwk] [out_prefix] [depth.tsv]
# Default taxonomy is SILVA (pass merged.tax.gtdb.csv as arg1 to use GTDB).
# Figures -> plots/ ; depth table in results/.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggtree)
  library(tidytree)
  library(treeio)
})

args      <- commandArgs(trailingOnly = TRUE)
tax_csv   <- ifelse(length(args) >= 1, args[1],
                    "results/sanger/merged/merged.tax.silva.csv")
tree_file <- ifelse(length(args) >= 2, args[2],
                    "results/sanger/merged/merged.reps.nwk")
out_pref  <- ifelse(length(args) >= 3, args[3],
                    "results/sanger/merged/merged.fasttree_genus_tree")
depth_tsv <- ifelse(length(args) >= 4, args[4],
                    "results/gourgouthakas_depth_table.tsv")

message("tax  : ", tax_csv, "  (taxonomy source)")
message("tree : ", tree_file, "  (de novo FastTree)")
message("out  : ", out_pref)

mode_chr <- function(x) {                 # most frequent non-empty value, or NA
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else names(which.max(table(x)))
}

# ---- 1. per-read GTDB genus/phylum -----------------------------------------
tax <- read_csv(tax_csv, show_col_types = FALSE) %>%
  mutate(genus = na_if(genus, "")) %>% filter(!is.na(genus))

genus_tot <- tax %>% group_by(genus) %>%
  summarise(n_isolates = n(), phylum = mode_chr(phylum), .groups = "drop")

# ---- 2. collapse the FastTree to one representative tip per genus -----------
tree <- ape::read.tree(tree_file)
tips <- tibble(label = tree$tip.label) %>%
  left_join(tax %>% select(query, genus, phylum, pct_id),
            by = c("label" = "query")) %>%
  filter(!is.na(genus))

# representative tip = the genus member with the highest GTDB % identity
reps <- tips %>% arrange(desc(pct_id)) %>% group_by(genus) %>%
  slice(1) %>% ungroup()
message("genera on the FastTree: ", nrow(reps),
        " (from ", length(tree$tip.label), " representative tips)")
if (nrow(reps) < 3) stop("Fewer than 3 genera; nothing to plot.")

pruned <- keep.tip(tree, reps$label)

tip_data <- tibble(label = pruned$tip.label) %>%
  left_join(reps %>% select(label, genus), by = "label") %>%
  left_join(genus_tot, by = "genus") %>%
  mutate(short_lab = genus, species = genus)

# ---- 3. depth table (placeholder; one row per displayed genus) --------------
# Non-destructive: keep existing rows (preserving real edits) and only append
# random placeholder rows for genera not yet present (e.g. SILVA genera missing
# from a GTDB-keyed table). Back up before changing.
depths <- c(0, -39, -220, -418, -678, -713, -900, -1050, -1100)
taxa   <- sort(unique(tip_data$short_lab))
rand_rows <- function(names) {
  set.seed(1)
  m <- matrix(sample(0:50, length(names) * length(depths), replace = TRUE),
              nrow = length(names), dimnames = list(NULL, as.character(depths)))
  bind_cols(tibble(taxon = names), as_tibble(m))
}
if (!file.exists(depth_tsv)) {
  dir.create(dirname(depth_tsv), recursive = TRUE, showWarnings = FALSE)
  write_tsv(rand_rows(taxa), depth_tsv)
  message("created placeholder depth table: ", depth_tsv)
} else {
  existing <- read_tsv(depth_tsv, show_col_types = FALSE)
  missing  <- setdiff(taxa, existing$taxon)
  if (length(missing) > 0) {
    file.copy(depth_tsv, paste0(depth_tsv, ".bak"), overwrite = TRUE)
    write_tsv(arrange(bind_rows(existing, rand_rows(missing)), taxon), depth_tsv)
    message("added ", length(missing), " missing genera to ", depth_tsv,
            " (backup -> ", basename(depth_tsv), ".bak)")
  }
}
depth_df <- read_tsv(depth_tsv, show_col_types = FALSE)

# ---- 4. tree (left) + total + depth table on one A4 page -------------------
fs_name <- 2.9; fs_cell <- 2.7; fs_hdr <- 3.4
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2",
                "#D55E00", "#CC79A7", "#000000", "#999999")

p <- ggtree(pruned) %<+% tip_data +
  geom_tippoint(aes(color = phylum, size = n_isolates))

dd <- p$data %>% filter(isTip)
xr <- max(dd$x); n <- nrow(dd)

x_name <- xr * 1.05
dx     <- xr * 0.32
x_tot  <- xr * 2.05
x0     <- x_tot + dx
hdr_y  <- n + 1

tot <- dd %>% distinct(short_lab, y, n_isolates) %>% mutate(xcell = x_tot)

long <- depth_df %>%
  pivot_longer(-taxon, names_to = "depth", values_to = "value") %>%
  mutate(depth = factor(depth, levels = as.character(depths)),
         xcell = x0 + (as.integer(depth) - 1) * dx) %>%
  inner_join(dd %>% select(short_lab, y), by = c("taxon" = "short_lab"))
long <- long %>%
  mutate(txtcol = if_else(value < 0.45 * max(value, na.rm = TRUE),
                          "white", "black"))
hdr <- tibble(depth = as.character(depths),
              xcell = x0 + (seq_along(depths) - 1) * dx)
x_right <- max(hdr$xcell) + dx

p_tab <- p +
  geom_text(data = dd, aes(x = x_name, y = y, label = short_lab, color = phylum),
            hjust = 0, size = fs_name, inherit.aes = FALSE) +
  geom_tile(data = tot, aes(x = xcell, y = y), fill = "grey92",
            width = dx * 0.92, height = 0.9, color = "grey80",
            linewidth = 0.1, inherit.aes = FALSE) +
  geom_text(data = tot, aes(x = xcell, y = y, label = n_isolates),
            size = fs_cell, fontface = "bold", inherit.aes = FALSE) +
  geom_tile(data = long, aes(x = xcell, y = y, fill = value),
            width = dx * 0.92, height = 0.9, color = "grey80",
            linewidth = 0.1, inherit.aes = FALSE) +
  geom_text(data = long, aes(x = xcell, y = y, label = value),
            color = long$txtcol, size = fs_cell, inherit.aes = FALSE) +
  geom_text(data = hdr, aes(x = xcell, y = hdr_y, label = depth),
            angle = 90, hjust = 0, size = fs_hdr, fontface = "bold",
            inherit.aes = FALSE) +
  annotate("text", x = x_name, y = hdr_y, label = "taxon", hjust = 0,
           fontface = "bold", size = fs_hdr) +
  annotate("text", x = x_tot, y = hdr_y, label = "total", angle = 90,
           hjust = 0, fontface = "bold", size = fs_hdr) +
  scale_fill_viridis_c(name = "abundance") +
  scale_color_manual(values = cb_palette, name = "phylum", na.value = "grey50") +
  scale_size_continuous(range = c(1, 5), name = "isolates") +
  coord_cartesian(xlim = c(0, x_right), ylim = c(0.5, n + 6), clip = "off") +
  theme_tree() +
  guides(
    fill  = guide_colorbar(title.position = "top", barheight = unit(0.3, "cm"),
                           barwidth = unit(3, "cm"), order = 1),
    color = guide_legend(title.position = "top", nrow = 2,
                         override.aes = list(size = 2.5), order = 2),
    size  = guide_legend(title.position = "top", nrow = 1, order = 3)) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "top",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.margin = margin(1, 1, 1, 1),
        plot.margin = margin(2, 6, 2, 4))

# ---- 5. write figures to plots/, data next to the input --------------------
dir.create("plots", showWarnings = FALSE)
fig_pref <- file.path("plots", basename(out_pref))
ggsave(paste0(fig_pref, ".pdf"), p_tab, device = "pdf",
       width = 21, height = 29.7, units = "cm", limitsize = FALSE)
ggsave(paste0(fig_pref, ".png"), p_tab, device = "png", dpi = 300,
       width = 21, height = 29.7, units = "cm", limitsize = FALSE)

p_circ <- ggtree(pruned, layout = "circular") %<+% tip_data +
  geom_tippoint(aes(color = genus, size = n_isolates)) +
  geom_tiplab(aes(label = short_lab), size = 2.2, align = TRUE,
              linesize = 0.1, offset = 0.02) +
  scale_color_viridis_d(name = "genus") +
  scale_size_continuous(range = c(1, 6)) +
  theme(legend.position = "right", plot.margin = margin(2, 2, 2, 2, "cm"))
ggsave(paste0(fig_pref, ".circular.png"), p_circ, device = "png", dpi = 300,
       height = 40, width = 45, units = "cm", limitsize = FALSE)

write.tree(pruned, file = paste0(out_pref, ".nwk"))
write_csv(tip_data, paste0(out_pref, ".tips.csv"))

message("done: figures -> ", fig_pref, ".{pdf,png,circular.png}")