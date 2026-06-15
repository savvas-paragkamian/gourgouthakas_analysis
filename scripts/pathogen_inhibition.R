#!/usr/bin/env Rscript

# ============================================================
# PART 1 - Publication figures for in vitro phytopathogen inhibition
# Input file: in_vitro_phytopathogens_inhibition.txt
#
# PART 2 (further below) - Ex-vivo Botrytis cinerea biocontrol
# statistics (per-d.p.i. ANOVA + Tukey HSD, AUDPC) and figures.
# Input file: ex-vivo-inhibition_B.c._SRL917_gourgouthakas.xlsx (sheet "raw")
# ============================================================

# -------------------------
# 1. Load packages
# -------------------------
packages <- c(
  "tidyverse",
  "readr",
  "readxl",
  "forcats",
  "viridis",
  "scales",
  "patchwork"
)

installed <- rownames(installed.packages())
for (p in packages) {
  if (!p %in% installed) install.packages(p)
}
lapply(packages, library, character.only = TRUE)

# -------------------------
# 2. Read data
# -------------------------
file_path <- "../data/in_vitro_phytopathogens_inhibition.txt"
file_taxonomy <- read_delim("../results/taxonomy_per_microbe.tsv", delim="\t") 

df <- read.delim(file_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
df$stab <- as.numeric(gsub("SRL","",df$Species))

df <- df |>
  left_join(file_taxonomy) |>
  mutate(gtdb_genus=if_else(is.na(gtdb_genus),Genus,gtdb_genus)) |>
  filter( is.na(silva_pct_id) | silva_pct_id>95)


# -------------------------
# 3. Define pathogen columns
# -------------------------
pathogen_cols <- c(
  "Xanthomonas campestris pv. Campestris_cm",
  "Paracidovorax citrulli_cm",
  "Ralstonia solanacearum_cm",
  "Clavibacter michiganensis_cm",
  "Verticillium dahliae_class",
  "Phytophthora nicotianae_class"
)

# -------------------------
# 4. Clean and standardize values
# -------------------------
df_clean <- df %>%
  mutate(across(all_of(pathogen_cols), ~case_when(
    . %in% c("NT", "__", "") ~ NA_character_,
    . == "Space-limiting" ~ "1",
    TRUE ~ as.character(.)
  ))) %>%
  mutate(across(all_of(pathogen_cols), as.numeric))

# -------------------------
# 5. Optional prettier pathogen labels
# -------------------------
pathogen_labels <- c(
  "Xanthomonas campestris pv. Campestris_cm" = "X. campestris",
  "Paracidovorax citrulli_cm"                = "P. citrulli",
  "Ralstonia solanacearum_cm"                = "R. solanacearum",
  "Clavibacter michiganensis_cm"             = "C. michiganensis",
  "Verticillium dahliae_class"               = "V. dahliae",
  "Phytophthora nicotianae_class"            = "P. nicotianae"
)

# -------------------------
# 6. Long format for plotting
# -------------------------
df_long <- df_clean %>%
  pivot_longer(
    cols = all_of(pathogen_cols),
    names_to = "Pathogen",
    values_to = "Inhibition"
  ) %>%
  mutate(
    Pathogen = recode(Pathogen, !!!pathogen_labels)
  )

# -------------------------
# 7. Publication-style theme
# -------------------------
theme_pub <- function(base_size = 12, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0),
      plot.subtitle = element_text(size = rel(0.95), hjust = 0),
      plot.caption = element_text(size = rel(0.85), colour = "grey30"),
      plot.margin = margin(10, 12, 10, 10)
    )
}

# -------------------------
# 8. Order isolates by total inhibition
# -------------------------
isolate_order <- df_long %>%
  group_by(Species) %>%
  summarise(total_signal = sum(Inhibition, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_signal)) %>%
  pull(Species)


# -------------------------
# 9. Heatmap
# -------------------------
heatmap_plot <- ggplot(df_long, aes(x = Pathogen, y = Species, fill = Inhibition)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_viridis(
    option = "magma",
    na.value = "grey92",
    name = "Inhibition\nscore",
    limits = c(0, 1.2),
    oob = squish
  ) +
  labs(
    title = "In vitro antagonistic activity of bacterial isolates",
    subtitle = "Heatmap of inhibition scores against bacterial and oomycete/fungal phytopathogens",
    x = NULL,
    y = "Isolate"
  ) +
  theme_pub(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 35, hjust = 1, vjust = 1
    ),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

# -------------------------
# 10. Summary for bar plot
#     Mean inhibition by genus across all tested pathogen-isolate values
# -------------------------
genus_summary <- df_long %>%
  group_by(gtdb_genus) %>%
  summarise(
    n = sum(!is.na(Inhibition)),
    mean_inhibition = mean(Inhibition, na.rm = TRUE),
    sd_inhibition = sd(Inhibition, na.rm = TRUE),
    se_inhibition = sd_inhibition / sqrt(n),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_inhibition)) %>%
  mutate(gtdb_genus = fct_reorder(gtdb_genus, mean_inhibition))

# -------------------------
# 11. Bar plot with error bars
#     Uses standard error; change to sd_inhibition if preferred
# -------------------------
bar_plot <- ggplot(genus_summary, aes(x = gtdb_genus, y = mean_inhibition)) +
  geom_col(width = 0.72, fill = "grey30") +
  geom_errorbar(
    aes(
      ymin = mean_inhibition - se_inhibition,
      ymax = mean_inhibition + se_inhibition
    ),
    width = 0.18,
    linewidth = 0.5
  ) +
  coord_flip() +
  labs(
    title = "Mean inhibition activity by genus",
    subtitle = "Bars show mean inhibition score; error bars represent standard error",
    x = NULL,
    y = "Mean inhibition score"
  ) +
  theme_pub(base_size = 12) + 
  theme(    axis.text.y = element_text(face="italic"))

# -------------------------
# 12. Save figures
# -------------------------
ggsave(
  filename = "../plots/figure_heatmap_inhibition.png",
  plot = heatmap_plot,
  width = 8.5,
  height = 11,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "../plots/figure_barplot_genus_inhibition.png",
  plot = bar_plot,
  width = 7,
  height = 4.8,
  dpi = 600,
  bg = "white"
)


############### Strong inhibutions per genus ###################


genus_summary <- df_long %>%
  group_by(gtdb_genus, Pathogen) %>%
  summarise(
    mean_inhibition = mean(Inhibition, na.rm = TRUE),
    strong_isolates = sum(Inhibition > 0.5, na.rm = TRUE),
    .groups = "drop"
  )

# Short pathogen labels
genus_summary$Pathogen <- recode(
  genus_summary$Pathogen,
  "Xanthomonas campestris pv. Campestris_cm"="Xanthomonas",
  "Paracidovorax citrulli_cm"="Paracidovorax",
  "Ralstonia solanacearum_cm"="Ralstonia",
  "Clavibacter michiganensis_cm"="Clavibacter",
  "Verticillium dahliae_class"="Verticillium",
  "Phytophthora nicotianae_class"="Phytophthora"
)

# Order genera
genus_summary <- genus_summary %>%
  mutate(Genus = fct_reorder(gtdb_genus, mean_inhibition, .fun = mean))

# Convert to factor for discrete colors
genus_summary$strong_isolates <- factor(genus_summary$strong_isolates)

# Color-blind safe discrete palette
cb_palette <- c(
"#56B4E9",
"#009E73",
"#E69F00",
"#D55E00",
"gray60",
"#CC79A7",
"#0072B2"
)

# -------------------------
# Plot
# -------------------------
p <- ggplot(
  genus_summary,
  aes(
    x = Pathogen,
    y = Genus,
    size = mean_inhibition,
    fill = strong_isolates
  )
) +
  geom_point(shape = 21, color="black", stroke=0.4) +
  scale_size(range = c(3,14), name="Mean inhibition") +
  scale_fill_manual(values = cb_palette, name="Isolates >0.5") +
  labs(
#    title = "Antagonistic activity of bacterial genera",
#    subtitle = "Bubble size = mean inhibition,\ncolor = number of strong inhibitory isolates (>0.5)",
    x = "Phytopathogen",
    y = "Bacterial genus"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1, face="italic"),
    axis.text.y = element_text(face="italic"),
    axis.title = element_text(face="bold"),
    legend.title = element_text(face="bold"),
    legend.text = element_text(face="italic"),
    plot.title = element_text(face="bold", size=14)
  )


  # -------------------------
# Save
# -------------------------
ggsave("../plots/bubble_strong_isolates.png", p, width=6, height=10, dpi=600)


# ============================================================
# ============================================================
# PART 2 - Ex-vivo Botrytis cinerea biocontrol statistics
#   Input: ../data/ex-vivo-inhibition_B.c._SRL917_gourgouthakas.xlsx
#          sheet "raw" (one row per replicate)
#   Columns: treatment | 1..6 d.p.i. disease severity (%) | AUDPC
#   Treatments: Control, B. cinerea, B.c.+X, B.c.+SRL917
#
#   Statistics:
#     - per-d.p.i. one-way ANOVA across treatments + Tukey HSD post-hoc
#     - one-way ANOVA + Tukey HSD on AUDPC
#   Figures:
#     - grouped (position-dodge) bar plot of mean disease severity per d.p.i.
#     - bar plot of mean AUDPC per treatment
#   All with SE error bars. Stats written to ../results/.
# ============================================================
# ============================================================

# -------------------------
# 2.1 Read the raw replicate-level sheet
# -------------------------
ex_path <- "../data/ex-vivo-inhibition_B.c._SRL917_gourgouthakas.xlsx"

ex_raw <- read_excel(ex_path, sheet = "raw")

# First (unnamed) column holds the treatment label
names(ex_raw)[1] <- "Treatment"

# d.p.i. severity columns, in day order
dpi_cols <- c("1 d.p.i.", "2 d.p.i.", "3 d.p.i.",
              "4 d.p.i.", "5 d.p.i.", "6 d.p.i.")

# Fixed, biologically meaningful treatment order
treatment_levels <- c("Control", "B. cinerea", "B.c.+X", "B.c.+SRL917")

ex_raw <- ex_raw %>%
  mutate(
    Treatment = factor(str_trim(Treatment), levels = treatment_levels),
    across(all_of(c(dpi_cols, "AUDPC")), as.numeric)
  ) %>%
  filter(!is.na(Treatment))

# -------------------------
# 2.1b AUDPC from the raw d.p.i. severities (trapezoidal rule)
#   Recomputes AUDPC so it does not depend on the sheet's precomputed column.
#     x = severity values for one replicate
#     t = matching time points (days); defaults to 1, 2, ... in column order
#   Trapezoidal integration: sum of (y_i + y_{i+1})/2 * (t_{i+1} - t_i).
#   NA-tolerant; needs >= 2 valid points or returns NA.
# -------------------------
audpc2 <- function(x, t = seq_along(x)) {
  ok <- !is.na(x) & !is.na(t)
  x <- x[ok]; t <- t[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(t); x <- x[o]; t <- t[o]
  sum((x[-1] + x[-length(x)]) / 2 * diff(t))
}

# day number parsed from each d.p.i. column name ("3 d.p.i." -> 3)
dpi_days <- readr::parse_number(dpi_cols)

# per-replicate AUDPC recomputed from the raw severities.
# The sheet's AUDPC anchors disease at inoculation (day 0, severity 0), i.e. it
# includes the leading day0->day1 trapezoid; we prepend (t=0, y=0) to match it.
ex_raw$AUDPC2 <- apply(as.matrix(ex_raw[dpi_cols]), 1,
                       function(r) audpc2(c(0, r), t = c(0, dpi_days)))

# sanity check vs the sheet's precomputed AUDPC column
cat("\n=== AUDPC recompute (audpc2) vs sheet AUDPC ===\n")
cat("max abs difference:",
    max(abs(ex_raw$AUDPC2 - ex_raw$AUDPC), na.rm = TRUE), "\n")

# -------------------------
# 2.2 Per-d.p.i. ANOVA + Tukey HSD
#     A column with no within-data variance (e.g. 1 d.p.i. is constant)
#     cannot be tested; it is reported as NA rather than crashing the run.
# -------------------------
anova_one <- function(dat, response) {
  vals <- dat[[response]]
  # need >1 group with data and some variance to run an ANOVA
  if (length(unique(na.omit(vals))) < 2 ||
      dplyr::n_distinct(dat$Treatment[!is.na(vals)]) < 2) {
    return(list(
      anova = tibble(dpi = response, df_between = NA_real_, df_within = NA_real_,
                     F = NA_real_, p = NA_real_,
                     note = "no variance - not testable"),
      tukey = NULL
    ))
  }
  fit  <- aov(reformulate("Treatment", response = response), data = dat)
  smry <- summary(fit)[[1]]
  aov_tbl <- tibble(
    dpi        = response,
    df_between = smry[["Df"]][1],
    df_within  = smry[["Df"]][2],
    F          = smry[["F value"]][1],
    p          = smry[["Pr(>F)"]][1],
    note       = NA_character_
  )
  tuk <- as.data.frame(TukeyHSD(fit)$Treatment)
  tuk <- tibble(
    dpi        = response,
    comparison = rownames(tuk),
    diff       = tuk$diff,
    lwr        = tuk$lwr,
    upr        = tuk$upr,
    p_adj      = tuk$`p adj`
  )
  list(anova = aov_tbl, tukey = tuk)
}

dpi_results <- lapply(dpi_cols, anova_one, dat = ex_raw)

anova_dpi <- bind_rows(lapply(dpi_results, `[[`, "anova"))
tukey_dpi <- bind_rows(lapply(dpi_results, `[[`, "tukey"))

# -------------------------
# 2.3 AUDPC ANOVA + Tukey HSD
# -------------------------
audpc_res   <- anova_one(ex_raw, "AUDPC")
anova_audpc <- audpc_res$anova %>% mutate(dpi = "AUDPC")
tukey_audpc <- audpc_res$tukey %>% mutate(dpi = "AUDPC")

# -------------------------
# 2.4 Write statistics to ../results/
# -------------------------
dir.create("../results", showWarnings = FALSE)

write_tsv(bind_rows(anova_dpi, anova_audpc),
          "../results/ex_vivo_anova.tsv")
write_tsv(bind_rows(tukey_dpi, tukey_audpc),
          "../results/ex_vivo_tukey.tsv")

cat("\n=== Ex-vivo per-d.p.i. ANOVA ===\n");  print(as.data.frame(anova_dpi))
cat("\n=== Ex-vivo AUDPC ANOVA ===\n");        print(as.data.frame(anova_audpc))
cat("\n=== Ex-vivo Tukey HSD (d.p.i.) ===\n"); print(as.data.frame(tukey_dpi))
cat("\n=== Ex-vivo Tukey HSD (AUDPC) ===\n");  print(as.data.frame(tukey_audpc))

# -------------------------
# 2.5 Summaries for plotting (mean +/- SE)
# -------------------------
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# Long over d.p.i. for the grouped bar plot
ex_long <- ex_raw %>%
  pivot_longer(all_of(dpi_cols), names_to = "dpi", values_to = "severity") %>%
  mutate(dpi = factor(dpi, levels = dpi_cols))

dpi_summary <- ex_long %>%
  group_by(Treatment, dpi) %>%
  summarise(
    n      = sum(!is.na(severity)),
    mean   = mean(severity, na.rm = TRUE),
    se     = se(severity),
    .groups = "drop"
  )

audpc_summary <- ex_raw %>%
  group_by(Treatment) %>%
  summarise(
    n      = sum(!is.na(AUDPC)),
    mean   = mean(AUDPC, na.rm = TRUE),
    se     = se(AUDPC),
    .groups = "drop"
  )

# Colour-blind safe palette (Okabe-Ito), one colour per treatment.
# Hues chosen for maximal separation under all colour-vision types:
# grey (neutral baseline) / red-orange (disease) / blue / bluish-green.
treatment_palette <- c(
  "Control"     = "#999999",  # grey
  "B. cinerea"  = "#D55E00",  # vermillion
  "B.c.+X"      = "#0072B2",  # blue
  "B.c.+SRL917" = "#009E73"   # bluish green
)

# -------------------------
# 2.6 Grouped (position-dodge) bar plot: disease severity per d.p.i.
#     Significance stars = each treatment vs B. cinerea (pathogen alone) at that
#     d.p.i. (Tukey HSD post-hoc p-adjusted): *** <0.001, ** <0.01, * <0.05, ns.
# -------------------------
dodge <- position_dodge(width = 0.8)

sig_stars <- function(p) cut(
  p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns"), right = FALSE
)

# Tukey comparisons that involve B. cinerea; the starred treatment is the other
# side of each pairwise contrast (no treatment label contains a "-", so the
# split is unambiguous). One row per (dpi, treatment).
ref_treatment <- "B. cinerea"
stars_vs_pathogen <- tukey_dpi %>%
  separate(comparison, into = c("side_a", "side_b"), sep = "-", remove = FALSE) %>%
  filter(side_a == ref_treatment | side_b == ref_treatment) %>%
  transmute(
    dpi       = dpi,
    Treatment = if_else(side_a == ref_treatment, side_b, side_a),
    stars     = as.character(sig_stars(p_adj))
  )

# Attach to the bar summary so each label sits just above its bar+SE.
# B. cinerea (the reference) keeps a blank label so position_dodge stays aligned.
stars_df <- dpi_summary %>%
  mutate(dpi = as.character(dpi)) %>%
  left_join(stars_vs_pathogen, by = c("dpi", "Treatment")) %>%
  mutate(
    stars     = ifelse(is.na(stars), "", stars),
    dpi       = factor(dpi, levels = dpi_cols),
    Treatment = factor(Treatment, levels = treatment_levels),
    y         = mean + se
  )

bar_dpi <- ggplot(dpi_summary, aes(x = dpi, y = mean, fill = Treatment)) +
  geom_col(width = 0.7, position = dodge, colour = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2, linewidth = 0.4, position = dodge
  ) +
  geom_text(
    data = stars_df,
    aes(x = dpi, y = y, label = stars, group = Treatment),
    position = dodge, vjust = -0.6, size = 3.2, fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = treatment_palette, name = "Treatment") +
  scale_y_continuous(
    breaks = scales::breaks_width(5),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Days post inoculation",
    y = "Spot diameter (mm)",
    caption = "Stars: treatment vs B. cinerea at each d.p.i. (Tukey HSD); *** p<0.001, ** p<0.01, * p<0.05, ns p>=0.05"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top",
        legend.text = element_text(face = "italic"))

# -------------------------
# 2.7 Bar plot: AUDPC per treatment
#     Stars = each treatment vs B. cinerea (Tukey HSD on AUDPC), above error bars.
# -------------------------
stars_audpc <- tukey_audpc %>%
  separate(comparison, into = c("side_a", "side_b"), sep = "-", remove = FALSE) %>%
  filter(side_a == ref_treatment | side_b == ref_treatment) %>%
  transmute(
    Treatment = if_else(side_a == ref_treatment, side_b, side_a),
    stars     = as.character(sig_stars(p_adj))
  ) %>%
  right_join(audpc_summary, by = "Treatment") %>%
  mutate(
    stars     = ifelse(is.na(stars), "", stars),  # B. cinerea reference = blank
    Treatment = factor(Treatment, levels = treatment_levels),
    y         = mean + se
  )

bar_audpc <- ggplot(audpc_summary,
                    aes(x = Treatment, y = mean, fill = Treatment)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.18, linewidth = 0.4
  ) +
  geom_text(
    data = stars_audpc,
    aes(x = Treatment, y = y, label = stars),
    vjust = -0.6, size = 3.5, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = treatment_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Ex vivo AUDPC by treatment",
    subtitle = "Area under the disease progress curve; error bars = SE",
    x = NULL,
    y = "AUDPC",
    caption = "vs B. cinerea (Tukey HSD): *** p<0.001, ** p<0.01, * p<0.05, ns"
  ) +
  theme_pub(base_size = 12) +
  theme(axis.text.x = element_text(face = "italic"))

# -------------------------
# 2.8 Save figures
# -------------------------
ggsave("../plots/ex_vivo_barplot_dpi.png",   bar_dpi,
       width = 8, height = 5, dpi = 600, bg = "white")
ggsave("../plots/ex_vivo_barplot_audpc.png", bar_audpc,
       width = 6, height = 5, dpi = 600, bg = "white")

