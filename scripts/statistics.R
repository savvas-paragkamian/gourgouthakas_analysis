#!/usr/bin/env Rscript

# ============================================================
# Publication figures for in vitro phytopathogen inhibition
# Input file: in_vitro_phytopathogens_inhibition.txt
# ============================================================

# -------------------------
# 1. Load packages
# -------------------------
packages <- c(
  "tidyverse",
  "readr",
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

df <- read.delim(file_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

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
  group_by(Genus) %>%
  summarise(
    n = sum(!is.na(Inhibition)),
    mean_inhibition = mean(Inhibition, na.rm = TRUE),
    sd_inhibition = sd(Inhibition, na.rm = TRUE),
    se_inhibition = sd_inhibition / sqrt(n),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_inhibition)) %>%
  mutate(Genus = fct_reorder(Genus, mean_inhibition))

# -------------------------
# 11. Bar plot with error bars
#     Uses standard error; change to sd_inhibition if preferred
# -------------------------
bar_plot <- ggplot(genus_summary, aes(x = Genus, y = mean_inhibition)) +
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
  group_by(Genus, Pathogen) %>%
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
  mutate(Genus = fct_reorder(Genus, mean_inhibition, .fun = mean))

# Convert to factor for discrete colors
genus_summary$strong_isolates <- factor(genus_summary$strong_isolates)

# Color-blind safe discrete palette
cb_palette <- c(
"#56B4E9",
"#009E73",
"#E69F00",
"#D55E00",
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
    title = "Antagonistic activity of bacterial genera",
    subtitle = "Bubble size = mean inhibition,\ncolor = number of strong inhibitory isolates (>0.5)",
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
ggsave("../plots/bubble_strong_isolates.png", p, width=6, height=8, dpi=600)

