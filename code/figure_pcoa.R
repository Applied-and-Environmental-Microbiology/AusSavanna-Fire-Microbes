
# Required packages
require(ggtext)
require(viridis)
require(tidyverse)

# (1) Setup --------------------------------------------------------------------

#### Data ####

# Define the factor levels
litter_levels <- c("Leaf", "Grass")
regime_levels <- c("E2", "L2", "U")
community_levels <- c("Soil bacteria", "Litter bacteria" ,"Soil fungi", "Litter fungi")

# Read in the sample metadata
metadata_fungi_soil <- data.table::fread("data/sample_metadata_fungi.txt") %>%
  filter(community == "Soil") %>%
  select(
    everything(), PCo1 = "pco1_54_soil_fungi", PCo2 = "pco2_39_soil_fungi",
    -c("community", "pco1_57_leaf_fungi", "pco2_38_leaf_fungi")
    ) %>%
  # Order the regime and litter type factors
  mutate(
    regime = factor(regime, levels = regime_levels),
    litter_type = factor(litter_type, levels = litter_levels)
  )
metadata_fungi_litter <- data.table::fread("data/sample_metadata_fungi.txt") %>%
  filter(community == "Leaf") %>%
  select(
    everything(), PCo1 = "pco1_57_leaf_fungi", PCo2 = "pco2_38_leaf_fungi",
    -c("community", "pco1_54_soil_fungi", "pco2_39_soil_fungi")
    ) %>%
  # Order the regime and litter type factors
  mutate(
    regime = factor(regime, levels = regime_levels),
    litter_type = factor(litter_type, levels = litter_levels)
  )
metadata_bacteria_soil <- data.table::fread("data/sample_metadata_bacteria.txt") %>%
  filter(community == "Soil") %>%
  select(
    everything(), PCo1 = "pco1_54_soil_bacteria", PCo2 = "pco2_40_soil_bacteria",
    -c("community", "pco1_49_leaf_bacteria", "pco2_42_leaf_bacteria")
    ) %>%
  # Order the regime and litter type factors
  mutate(
    regime = factor(regime, levels = regime_levels),
    litter_type = factor(litter_type, levels = litter_levels)
  )
metadata_bacteria_litter <- data.table::fread("data/sample_metadata_bacteria.txt") %>%
  filter(community == "Leaf") %>%
  select(
    everything(), PCo1 = "pco1_49_leaf_bacteria", PCo2 = "pco2_42_leaf_bacteria",
    -c("community", "pco1_54_soil_bacteria", "pco2_40_soil_bacteria")
    ) %>%
  # Order the regime and litter type factors
  mutate(
    regime = factor(regime, levels = regime_levels),
    litter_type = factor(litter_type, levels = litter_levels)
  )

# Format the axis labels
axis_labels <- c(
  "pco1_54_soil_bacteria" = "PCo1 (54%)",
  "pco2_40_soil_bacteria" = "PCo2 (40%)",
  "pco1_49_leaf_bacteria" = "PCo1 (49%)",
  "pco2_42_leaf_bacteria" = "PCo2 (42%)",
  "pco1_57_leaf_fungi" = "PCo1 (57%)",
  "pco2_38_leaf_fungi" = "PCo2 (38%)",
  "pco1_54_soil_fungi" = "PCo1 (54%)",
  "pco2_39_soil_fungi" = "PCo2 (39%)"
)

#### Theme ####

tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 8

my_theme <- theme_minimal() +
  theme(
    panel.border = element_rect(linewidth = 0.5, colour = "grey70", fill = NA),
    panel.grid = element_blank(),
    axis.title = element_text(size = title_size),
    axis.text = element_text(size = text_size),
    plot.title = element_text(hjust = 0.5, size = strip_size, face = "bold"),
    plot.tag = element_markdown(size = tag_size, face = "bold"),
    legend.position = "none",
    aspect.ratio = 1
  )

#### Dummy legend ####

# Define colours
my_cols <- c("E2" = "#3B0F70CC", "L2" = "#B63679CC", "U" = "#FE9F6DCC")


# Dummy tibble with all regime/shape combinations
dummy_data <- expand.grid(
  regime = names(my_cols),
  litter_type = c("Leaf", "Grass")
)

# Dummy plot that only serves to generate legend
dummy_plot <- ggplot(dummy_data, aes(x = 1, y = 1, colour = regime, shape = litter_type)) +
  geom_point(size = 3) +
  scale_colour_manual(values = my_cols) +
  scale_shape_manual(values = c("Leaf" = 16, "Grass" = 1)) +
  labs(colour = "Regime", shape = "Litter type") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = title_size),
    legend.text = element_text(size = text_size)
  )

# Extract the legend
legend_dummy <- cowplot::get_legend(dummy_plot)

# (2) Crate plots --------------------------------------------------------------

# Create the figure
figure <- patchwork::wrap_plots(
  ggplot(
    data = metadata_bacteria_soil,
    aes(x = PCo1, y = PCo2, colour = regime, fill = regime, shape = litter_type)
  ) +
    geom_point(size = 2.2, alpha = 0.7) +
    # Elipse for leaf litter
    stat_ellipse(
      data = metadata_bacteria_soil %>%
        filter(litter_type == "Leaf"),
      aes(x = PCo1, y = PCo2, group = regime, colour = regime),
      level = 0.95, linewidth = 0.3, alpha = 0.6,
      show.legend = FALSE
    ) +
    # Elipse for grass litter
    stat_ellipse(
      data = metadata_bacteria_soil %>%
        filter(litter_type == "Grass"),
      aes(x = PCo1, y = PCo2, group = regime, colour = regime),
      level = 0.95, linewidth = 0.3, alpha = 0.6, linetype = "dashed",
      show.legend = FALSE
    ) +
    scale_colour_manual(values = my_cols) +
    scale_fill_manual(values = my_cols) +
    scale_shape_manual(values = c("Leaf" = 16, "Grass" = 1)) +
    labs(
      x = axis_labels["pco1_54_soil_bacteria"],
      y = axis_labels["pco2_40_soil_bacteria"],
      colour = "Regime",
      shape = "Litter type",
      title = "Soil",
      tag = "(**a**)"
    ) +
    my_theme +
   theme(
     # Adjust the position of the tag to account for the strip
     plot.tag.position = c(0.04, 0.95),
     # Adjust the margin
     plot.margin = margin(t = 0, r = 0, b = 5, l = 0),
     ),
  
  ggplot(
    data = metadata_bacteria_litter,
    aes(x = PCo1, y = PCo2, colour = regime, fill = regime, shape = litter_type)
  ) +
    geom_point(size = 2.2, alpha = 0.7) +
    # Elipse for leaf litter
    stat_ellipse(
      data = metadata_bacteria_litter %>%
        filter(litter_type == "Leaf"),
      aes(x = PCo1, y = PCo2, group = regime, colour = regime),
      level = 0.95, linewidth = 0.3, alpha = 0.6,
      show.legend = FALSE
    ) +
    # Elipse for grass litter
    stat_ellipse(
      data = metadata_bacteria_litter %>%
        filter(litter_type == "Grass"),
      aes(x = PCo1, y = PCo2, group = regime, colour = regime),
      level = 0.95, linewidth = 0.3, alpha = 0.6, linetype = "dashed",
      show.legend = FALSE
    ) +
    scale_colour_manual(values = my_cols) +
    scale_fill_manual(values = my_cols) +
    scale_shape_manual(values = c("Leaf" = 16, "Grass" = 1)) +
    labs(
      x = axis_labels["pco1_49_leaf_bacteria"],
      y = axis_labels["pco2_42_leaf_bacteria"],
      colour = "Regime",
      shape = "Litter type",
      title = "Litter"
    )  +
    my_theme +
    theme(
      # Adjust the margin
      plot.margin = margin(t = 0, r = 0, b = 5, l = 5),
    ),
  
  ggplot(
    data = metadata_fungi_soil,
    aes(x = PCo1, y = PCo2, colour = regime, fill = regime, shape = litter_type)
  ) +
    geom_point(size = 2.2, alpha = 0.7) +
    # Elipse for leaf litter
    stat_ellipse(
      data = metadata_fungi_soil %>%
        filter(litter_type == "Leaf"),
      aes(x = PCo1, y = PCo2, group = regime, colour = regime),
      level = 0.95, linewidth = 0.3, alpha = 0.6,
      show.legend = FALSE
    ) +
    # Elipse for grass litter
    stat_ellipse(
      data = metadata_fungi_soil %>%
        filter(litter_type == "Grass"),
      aes(x = PCo1, y = PCo2, group = regime, colour = regime),
      level = 0.95, linewidth = 0.3, alpha = 0.6, linetype = "dashed",
      show.legend = FALSE
    ) +
    scale_colour_manual(values = my_cols) +
    scale_fill_manual(values = my_cols) +
    scale_shape_manual(values = c("Leaf" = 16, "Grass" = 1)) +
    labs(
      x = axis_labels["pco1_54_soil_fungi"],
      y = axis_labels["pco2_39_soil_fungi"],
      colour = "Regime",
      shape = "Litter type",
      tag = "(**b**)"
    ) +
    my_theme +
    theme(
      # Adjust the margin
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    ),
  
  ggplot(
    data = metadata_fungi_litter,
    aes(x = PCo1, y = PCo2, colour = regime, fill = regime, shape = litter_type)
  ) +
    geom_point(size = 2.2, alpha = 0.7) +
    # Elipse for leaf litter
    stat_ellipse(
      data = metadata_fungi_litter %>%
        filter(litter_type == "Leaf"),
      aes(x = PCo1, y = PCo2, group = regime, colour = regime),
      level = 0.95, linewidth = 0.2, alpha = 0.6,
      show.legend = FALSE
    ) +
    # Elipse for grass litter
    stat_ellipse(
      data = metadata_fungi_litter %>%
        filter(litter_type == "Grass"),
      aes(x = PCo1, y = PCo2, group = regime, colour = regime),
      level = 0.95, linewidth = 0.3, alpha = 0.6, linetype = "dashed",
      show.legend = FALSE
    ) +
    scale_colour_manual(values = my_cols) +
    scale_fill_manual(values = my_cols) +
    scale_shape_manual(values = c("Leaf" = 16, "Grass" = 1)) +
    labs(
      x = axis_labels["pco1_57_leaf_fungi"],
      y = axis_labels["pco2_38_leaf_fungi"],
      colour = "Regime",
      shape = "Litter type"
    ) +
    my_theme +
    theme(
      # Adjust the margin
      plot.margin = margin(t = 0, r = 0, b = 0, l = 5),
    )
)

# Add the legend
figure_final <- cowplot::plot_grid(
  figure, legend_dummy, rel_widths = c(1, 0.15)
)

# Save figure
ggsave(
  "output/figure_pcoa.pdf",
  width = 15, height = 13, units = "cm",
)

