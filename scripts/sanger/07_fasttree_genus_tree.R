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
# match the per-taxonomy depth table (08 writes .gtdb.tsv / .silva.tsv) to the
# taxonomy source so genus counts match the tip labels.
tax_db    <- ifelse(grepl("silva", tax_csv), "silva", "gtdb")
depth_tsv <- ifelse(length(args) >= 4, args[4],
                    paste0("results/gourgouthakas_depth_table.", tax_db, ".tsv"))

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

# depth table (built by 08): keep only genera with at least one isolate in it
# (0-filled where a genus is absent at a depth); drop the rest, no fabrication.
depths <- c(0, -39, -220, -418, -678, -713, -900, -1050, -1100)
if (!file.exists(depth_tsv))
  stop("depth table not found: ", depth_tsv, " -- run 08_taxonomy_table.py first")
depth_df <- read_tsv(depth_tsv, show_col_types = FALSE)
dropped <- setdiff(reps$genus, depth_df$taxon)
reps    <- reps %>% filter(genus %in% depth_df$taxon)
if (length(dropped) > 0)
  message("dropped ", length(dropped), " genera absent from ",
          basename(depth_tsv), " (no depth metadata): ",
          paste(head(dropped, 6), collapse = ", "),
          if (length(dropped) > 6) ", ..." else "")

message("genera on the FastTree: ", nrow(reps),
        " (from ", length(tree$tip.label), " representative tips)")
if (nrow(reps) < 3) stop("Fewer than 3 genera; nothing to plot.")

pruned <- keep.tip(tree, reps$label)

tip_data <- tibble(label = pruned$tip.label) %>%
  left_join(reps %>% select(label, genus), by = "label") %>%
  left_join(genus_tot, by = "genus") %>%
  mutate(short_lab = genus, species = genus)

# Count only the supplied-table isolates: the "total" column and the tip size
# are the depth-table row sum, so total == sum of the depth cells (inner join).
tip_data <- tip_data %>% select(-n_isolates) %>%
  left_join(depth_df %>%
              transmute(genus = taxon,
                        n_isolates = rowSums(across(all_of(as.character(depths))))),
            by = "genus")

# ---- 4. tree (left) + total + depth table on one A4 page -------------------
fs_name <- 2.9; fs_cell <- 2.7; fs_hdr <- 3.4
# Okabe-Ito palette ordered for maximum divergence between the few phyla shown,
# excluding the blue (#0072B2) used by the abundance heatmap below.
cb_palette <- c("#D55E00", "#009E73", "#E69F00", "#CC79A7",
                "#56B4E9", "#000000", "#999999")

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
  scale_fill_gradient(low = "white", high = "#0072B2", name = "abundance") +
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
ggsave(paste0(fig_pref, ".pdf"), p_tab, device = cairo_pdf,
       width = 21, height = 29.7, units = "cm", limitsize = FALSE)
ggsave(paste0(fig_pref, ".png"), p_tab, device = "png", dpi = 300,
       width = 21, height = 29.7, units = "cm", limitsize = FALSE)

p_circ <- ggtree(pruned, layout = "circular") %<+% tip_data +
  geom_tippoint(aes(color = genus, size = n_isolates)) +
  geom_tiplab(aes(label = short_lab), size = 2.2, align = TRUE,
              linesize = 0.1, offset = 0.02, fontface = "italic",
              family = "Ubuntu") +
  scale_color_viridis_d(name = "genus") +
  scale_size_continuous(range = c(1, 6)) +
  theme(legend.position = "right", plot.margin = margin(2, 2, 2, 2, "cm"))
ggsave(paste0(fig_pref, ".circular.png"), p_circ, device = "png", dpi = 300,
       height = 40, width = 45, units = "cm", limitsize = FALSE)

write.tree(pruned, file = paste0(out_pref, ".nwk"))
write_csv(tip_data, paste0(out_pref, ".tips.csv"))

message("done: figures -> ", fig_pref, ".{pdf,png,circular.png}")