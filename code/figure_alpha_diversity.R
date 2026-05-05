# Read in necessary libraries
require(ggtext)
require(tidyverse)

# Read in the data
load("data/alpha_diversity_plot_data.RData")

#### Theme ####
bacteria_colours <- c("Broadleaf" = "#4a7aab", "Grass" = "#b0c4d8")
fungi_colours    <- c("Broadleaf" = "#c4622d", "Grass" = "#e8c4a8")
organism_litter_colours <- c(
  "Grass (bacteria)"     = "#b0c4d8",
  "Broadleaf (bacteria)" = "#4a7aab",
  "Grass (fungi)"        = "#e8c4a8",
  "Broadleaf (fungi)"    = "#c4622d"
)
tag_size <- 14
strip_size <- 12
title_size <- 10
text_size <- 8
my_theme <- theme_minimal() +
  theme(
    panel.border = element_rect(linewidth = 0.5, colour = "grey70", fill = NA),
    panel.grid = element_line(colour = "grey90", linewidth = 0.15),
    axis.title = element_text(size = title_size),
    axis.text = element_text(size = text_size),
    strip.text = element_text(size = strip_size, face = "bold"),
    plot.tag = element_markdown(size = tag_size),
    legend.position = "none",
    aspect.ratio = 1
  )

#### Create dummy legends ####
dummy_bacteria_data <- data.frame(
  litter_type = factor(c("Broadleaf", "Grass"), levels = c("Broadleaf", "Grass")),
  value = 1:2
)
dummy_fungi_data <- data.frame(
  litter_type = factor(c("Broadleaf", "Grass"), levels = c("Broadleaf", "Grass")),
  value = 1:2
)
dummy_combined_data <- data.frame(
  group = factor(
    c("Grass (bacteria)", "Broadleaf (bacteria)", "Grass (fungi)", "Broadleaf (fungi)"),
    levels = c("Grass (bacteria)", "Broadleaf (bacteria)", "Grass (fungi)", "Broadleaf (fungi)")
  ),
  value = 1:4
)

dummy_legend_bacteria <- cowplot::get_legend(
  ggplot(dummy_bacteria_data, aes(x = value, y = value, fill = litter_type)) +
    geom_point(colour = alpha('grey30', 0.5), shape = 21, size = 5) +
    scale_fill_manual(values = bacteria_colours, name = "Litter type") +
    theme_void() +
    scale_y_continuous(limits = c(0, 0))
)

dummy_legend_fungi <- cowplot::get_legend(
  ggplot(dummy_fungi_data, aes(x = value, y = value, fill = litter_type)) +
    geom_point(colour = alpha('grey30', 0.5), shape = 21, size = 5) +
    scale_fill_manual(values = fungi_colours, name = "Litter type") +
    theme_void() +
    scale_y_continuous(limits = c(0, 0))
)

dummy_legend_combined <- cowplot::get_legend(
  ggplot(dummy_combined_data, aes(x = value, y = value, fill = group)) +
    geom_point(colour = alpha('grey30', 0.5), shape = 21, size = 5) +
    scale_fill_manual(values = organism_litter_colours, name = "Litter type") +
    theme_void() +
    scale_y_continuous(limits = c(0, 0))
)

#### Helper function for formatting intercepts ####
format_intercept <- function(mean_val, lower_val, upper_val) {
  format_val <- function(x) {
    if (abs(x) >= 100) {
      sprintf("%.0f", x)
    } else if (abs(x) >= 10) {
      sprintf("%.1f", x)
    } else {
      sprintf("%.2f", x)
    }
  }

  paste0(
    "&beta;<sub>0</sub> = ",
    format_val(mean_val),
    " [",
    format_val(lower_val),
    ", ",
    format_val(upper_val),
    "]"
  )
}

#### Prepare intercept annotations ####
intercept_annotations <- boot_intercepts %>%
  separate(Intercept, into = c("mean", "ci"), sep = " \\[", remove = FALSE) %>%
  separate(ci, into = c("lower", "upper"), sep = ", ", remove = FALSE) %>%
  mutate(
    mean = as.numeric(mean),
    lower = as.numeric(lower),
    upper = as.numeric(gsub("\\]", "", upper)),
    label = pmap_chr(list(mean, lower, upper), format_intercept),
    x = -Inf,
    y = Inf,
    hjust = 1,
    vjust = 0.1
  )

#### FIGURE 1: Shannon diversity (Bacteria & Fungi) ####

## Bacteria Shannon diversity plot
bacteria_shannon_plot <- ggplot() +
  geom_point(
    data = raw_data %>%
      filter(Response %in% c("Soil bacteria shannon", "Litter bacteria shannon")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria shannon", "Litter bacteria shannon"))),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.5),
    alpha = 0.5
  ) +
  geom_errorbar(
    data = emmean %>%
      filter(Response %in% c("Soil bacteria shannon", "Litter bacteria shannon")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria shannon", "Litter bacteria shannon"))),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>%
      filter(Response %in% c("Soil bacteria shannon", "Litter bacteria shannon")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria shannon", "Litter bacteria shannon"))),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    plot.tag.position = c(0.03, 0.92)
  ) +
  scale_y_continuous(
    limits = c(4, 6.5),
    breaks = scales::pretty_breaks(n = 3)
  ) +
  scale_color_manual(values = bacteria_colours) +
  scale_fill_manual(values = bacteria_colours) +
  labs(
    x = NULL,
    y = "Shannon diversity",
    tag = "(**a**)"
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil bacteria shannon" = "Soil",
      "Litter bacteria shannon" = "Litter"
    ))
  )

## Bacteria Shannon effect plot
bacteria_shannon_density_nudged <- boot_densities %>%
  filter(Response %in% c("Soil bacteria shannon", "Litter bacteria shannon")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil bacteria shannon", "Litter bacteria shannon")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

bacteria_shannon_effect_nudged <- boot_effects %>%
  filter(Response %in% c("Soil bacteria shannon", "Litter bacteria shannon")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil bacteria shannon", "Litter bacteria shannon")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

bacteria_shannon_effect_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggridges::geom_density_ridges(
    data = bacteria_shannon_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type,
        group = interaction(regime, litter_type)),
    stat = "identity",
    scale = 0.75,
    colour = alpha("grey30", 0.5),
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = bacteria_shannon_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = bacteria_shannon_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  geom_richtext(
    data = intercept_annotations %>%
      filter(Response %in% c("Soil bacteria shannon", "Litter bacteria shannon")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria shannon", "Litter bacteria shannon"))),
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    size = 3,
    fill = NA,
    label.color = NA
  ) +
  scale_y_continuous(
    breaks = unique(bacteria_shannon_density_nudged$regime_numeric),
    labels = unique(bacteria_shannon_density_nudged$regime),
    limits = c(0.5, max(bacteria_shannon_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(values = bacteria_colours) +
  coord_flip() +
  my_theme +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(
    x = "Effect size",
    y = NULL
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil bacteria shannon" = "Soil",
      "Litter bacteria shannon" = "Litter"
    ))
  )

## Fungi Shannon diversity plot
fungi_shannon_plot <- ggplot() +
  geom_point(
    data = raw_data %>%
      filter(Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi shannon", "Litter fungi shannon"))),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.5),
    alpha = 0.5
  ) +
  geom_errorbar(
    data = emmean %>%
      filter(Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi shannon", "Litter fungi shannon"))),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>%
      filter(Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi shannon", "Litter fungi shannon"))),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank()
  ) +
  scale_y_continuous(
    limits = c(1.5, 5),
    breaks = scales::pretty_breaks(n = 3)
  ) +
  scale_color_manual(values = fungi_colours) +
  scale_fill_manual(values = fungi_colours) +
  labs(
    x = NULL,
    y = "Shannon diversity",
    tag = "(**b**)"
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil fungi shannon" = "Soil",
      "Litter fungi shannon" = "Litter"
    ))
  )

## Fungi Shannon effect plot
fungi_shannon_density_nudged <- boot_densities %>%
  filter(Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil fungi shannon", "Litter fungi shannon")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

fungi_shannon_effect_nudged <- boot_effects %>%
  filter(Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil fungi shannon", "Litter fungi shannon")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

fungi_shannon_effect_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggridges::geom_density_ridges(
    data = fungi_shannon_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type,
        group = interaction(regime, litter_type)),
    stat = "identity",
    scale = 0.75,
    colour = alpha("grey30", 0.5),
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = fungi_shannon_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = fungi_shannon_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  geom_richtext(
    data = intercept_annotations %>%
      filter(Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi shannon", "Litter fungi shannon"))),
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    size = 3,
    fill = NA,
    label.color = NA
  ) +
  scale_y_continuous(
    breaks = unique(fungi_shannon_density_nudged$regime_numeric),
    labels = unique(fungi_shannon_density_nudged$regime),
    limits = c(0.5, max(fungi_shannon_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(values = fungi_colours) +
  coord_flip() +
  my_theme +
  theme(strip.text = element_blank()) +
  labs(
    x = "Effect size",
    y = NULL
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil fungi shannon" = "Soil",
      "Litter fungi shannon" = "Litter"
    ))
  )

## Assemble Figure 1
figure_shannon_diversity <- patchwork::wrap_plots(
  bacteria_shannon_plot, bacteria_shannon_effect_plot,
  fungi_shannon_plot, fungi_shannon_effect_plot,
  ncol = 1
)

figure_shannon_diversity_final <- cowplot::plot_grid(
  figure_shannon_diversity, dummy_legend_combined, rel_widths = c(1, 0.24)
)

figure_shannon_raw_only <- patchwork::wrap_plots(
  bacteria_shannon_plot +
    theme(
      legend.position = "right",
      legend.justification = c(0, 0.5), 
      plot.tag.position = c(0.03, 0.92) 
    ) +
    guides(
      fill   = guide_legend(title = "Bacteria litter type"),
      colour = guide_legend(title = "Bacteria litter type")
    ),
  fungi_shannon_plot +
    theme(
      legend.position = "right",
      legend.justification = c(0, 0.5),
      axis.text.x = element_text(size = text_size)
    ) +
    guides(
      fill   = guide_legend(title = "Fungi litter type"),
      colour = guide_legend(title = "Fungi litter type")
    ),
  ncol = 1
)

figure_shannon_effect_only <- patchwork::wrap_plots(
  bacteria_shannon_effect_plot +
    labs(tag = "(**a**)") +
    theme(
      legend.position = "right",
      legend.justification = c(0, 0.5),
      strip.text = element_text(size = strip_size, face = "bold", vjust = 1),
      plot.tag.position = c(0.03, 0.97)
    ) +
    guides(fill = guide_legend(title = "Bacteria litter type")),
  fungi_shannon_effect_plot +
    labs(tag = "(**b**)") +
    theme(
      legend.position = "right",
      legend.justification = c(0, 0.5)
    ) +
    guides(fill = guide_legend(title = "Fungi litter type")),
  ncol = 1
)

# Create output directory if it doesn't exist
if (!dir.exists("output")) {
  dir.create("output")
}

# Save Figure 1 — combined
ggsave(
  filename = "output/figure_shannon_diversity.tif",
  plot = figure_shannon_diversity_final,
  width = 13.5,
  height = 20,
  units = "cm",
  bg = "white"
)
ggsave(
  filename = "output/figure_shannon_diversity.png",
  plot = figure_shannon_diversity_final,
  width = 13.5,
  height = 20,
  units = "cm",
  bg = "white",
  dpi = 300
)

# Save Figure 1 — raw values and marginal means only
ggsave(
  filename = "output/figure_shannon_diversity_raw_only.tif",
  plot = figure_shannon_raw_only,
  width = 13.5,
  height = 9,
  units = "cm",
  bg = "white"
)
ggsave(
  filename = "output/figure_shannon_diversity_raw_only.png",
  plot = figure_shannon_raw_only,
  width = 13.5,
  height = 9,
  units = "cm",
  bg = "white",
  dpi = 300
)

# Save Figure 1 — effect sizes only
ggsave(
  filename = "output/figure_shannon_diversity_effect_only.tif",
  plot = figure_shannon_effect_only,
  width = 13.5,
  height = 9,
  units = "cm",
  bg = "white"
)
ggsave(
  filename = "output/figure_shannon_diversity_effect_only.png",
  plot = figure_shannon_effect_only,
  width = 13.5,
  height = 9,
  units = "cm",
  bg = "white",
  dpi = 300
)

#### FIGURE 2: Bacteria richness & evenness ####

## Bacteria richness plot
bacteria_richness_plot <- ggplot() +
  geom_point(
    data = raw_data %>%
      filter(Response %in% c("Soil bacteria richness", "Litter bacteria richness")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria richness", "Litter bacteria richness"))),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.5),
    alpha = 0.5
  ) +
  geom_errorbar(
    data = emmean %>%
      filter(Response %in% c("Soil bacteria richness", "Litter bacteria richness")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria richness", "Litter bacteria richness"))),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>%
      filter(Response %in% c("Soil bacteria richness", "Litter bacteria richness")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria richness", "Litter bacteria richness"))),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    plot.tag.position = c(0.03, 0.92)
  ) +
  scale_y_continuous(
    limits = c(
      0,
      raw_data %>%
        filter(Response %in% c("Soil bacteria richness", "Litter bacteria richness")) %>%
        pull(observed_value) %>%
        max()
    ),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = bacteria_colours) +
  scale_fill_manual(values = bacteria_colours) +
  labs(
    x = NULL,
    y = "Richness",
    tag = "(**a**)"
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil bacteria richness" = "Soil",
      "Litter bacteria richness" = "Litter"
    ))
  )

## Bacteria richness effect plot
bacteria_richness_density_nudged <- boot_densities %>%
  filter(Response %in% c("Soil bacteria richness", "Litter bacteria richness")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil bacteria richness", "Litter bacteria richness")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

bacteria_richness_effect_nudged <- boot_effects %>%
  filter(Response %in% c("Soil bacteria richness", "Litter bacteria richness")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil bacteria richness", "Litter bacteria richness")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

bacteria_richness_effect_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggridges::geom_density_ridges(
    data = bacteria_richness_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type,
        group = interaction(regime, litter_type)),
    stat = "identity",
    scale = 0.75,
    colour = alpha("grey30", 0.5),
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = bacteria_richness_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = bacteria_richness_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  geom_richtext(
    data = intercept_annotations %>%
      filter(Response %in% c("Soil bacteria richness", "Litter bacteria richness")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria richness", "Litter bacteria richness"))),
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    size = 3,
    fill = NA,
    label.color = NA
  ) +
  scale_y_continuous(
    breaks = unique(bacteria_richness_density_nudged$regime_numeric),
    labels = unique(bacteria_richness_density_nudged$regime),
    limits = c(0.5, max(bacteria_richness_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(values = bacteria_colours) +
  coord_flip() +
  my_theme +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(
    x = "Effect size",
    y = NULL
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil bacteria richness" = "Soil",
      "Litter bacteria richness" = "Litter"
    ))
  )

## Bacteria evenness plot
bacteria_evenness_plot <- ggplot() +
  geom_point(
    data = raw_data %>%
      filter(Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria evenness", "Litter bacteria evenness"))),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.5),
    alpha = 0.5
  ) +
  geom_errorbar(
    data = emmean %>%
      filter(Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria evenness", "Litter bacteria evenness"))),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>%
      filter(Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria evenness", "Litter bacteria evenness"))),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0.75, 1),
    breaks = scales::pretty_breaks(n = 3)
  ) +
  scale_color_manual(values = bacteria_colours) +
  scale_fill_manual(values = bacteria_colours) +
  labs(
    x = NULL,
    y = "Evenness",
    tag = "(**b**)"
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil bacteria evenness" = "Soil",
      "Litter bacteria evenness" = "Litter"
    ))
  )

## Bacteria evenness effect plot
bacteria_evenness_density_nudged <- boot_densities %>%
  filter(Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil bacteria evenness", "Litter bacteria evenness")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

bacteria_evenness_effect_nudged <- boot_effects %>%
  filter(Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil bacteria evenness", "Litter bacteria evenness")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

bacteria_evenness_effect_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggridges::geom_density_ridges(
    data = bacteria_evenness_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type,
        group = interaction(regime, litter_type)),
    stat = "identity",
    scale = 0.85,
    colour = alpha("grey30", 0.5),
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = bacteria_evenness_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = bacteria_evenness_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  geom_richtext(
    data = intercept_annotations %>%
      filter(Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
      mutate(Response = factor(Response, levels = c("Soil bacteria evenness", "Litter bacteria evenness"))),
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    size = 3,
    fill = NA,
    label.color = NA
  ) +
  scale_y_continuous(
    breaks = unique(bacteria_evenness_density_nudged$regime_numeric),
    labels = unique(bacteria_evenness_density_nudged$regime),
    limits = c(0.5, max(bacteria_evenness_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(values = bacteria_colours) +
  coord_flip() +
  my_theme +
  theme(strip.text = element_blank()) +
  labs(
    x = "Effect size",
    y = NULL
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil bacteria evenness" = "Soil",
      "Litter bacteria evenness" = "Litter"
    ))
  )

## Assemble Figure 2
figure_bacteria_richness_evenness <- patchwork::wrap_plots(
  bacteria_richness_plot, bacteria_richness_effect_plot,
  bacteria_evenness_plot, bacteria_evenness_effect_plot,
  ncol = 1
)

figure_bacteria_richness_evenness_final <- cowplot::plot_grid(
  figure_bacteria_richness_evenness, dummy_legend_bacteria, rel_widths = c(1, 0.24)
)

figure_bacteria_richness_evenness_raw_only <- cowplot::plot_grid(
  patchwork::wrap_plots(
    bacteria_richness_plot,
    bacteria_evenness_plot +
      theme(axis.text.x = element_text(size = text_size)),
    ncol = 1
  ),
  dummy_legend_bacteria,
  rel_widths = c(1, 0.24)
)

figure_bacteria_richness_evenness_effect_only <- cowplot::plot_grid(
  patchwork::wrap_plots(
    bacteria_richness_effect_plot +
      labs(tag = "(**a**)") +
      theme(
        strip.text = element_text(size = strip_size, face = "bold", vjust = 1),
        plot.tag.position = c(0.03, 0.97)
      ),
    bacteria_evenness_effect_plot +
      labs(tag = "(**b**)"),
    ncol = 1
  ),
  dummy_legend_bacteria,
  rel_widths = c(1, 0.24)
)

# Save Figure 2 — combined
ggsave(
  filename = "output/figure_bacteria_richness_evenness.tif",
  plot = figure_bacteria_richness_evenness_final,
  width = 13.5,
  height = 20,
  units = "cm",
  bg = "white"
)
ggsave(
  filename = "output/figure_bacteria_richness_evenness.png",
  plot = figure_bacteria_richness_evenness_final,
  width = 13.5,
  height = 20,
  units = "cm",
  bg = "white",
  dpi = 300
)

# Save Figure 2 — raw values and marginal means only
ggsave(
  filename = "output/figure_bacteria_richness_evenness_raw_only.tif",
  plot = figure_bacteria_richness_evenness_raw_only,
  width = 13.5,
  height = 10,
  units = "cm",
  bg = "white"
)
ggsave(
  filename = "output/figure_bacteria_richness_evenness_raw_only.png",
  plot = figure_bacteria_richness_evenness_raw_only,
  width = 13.5,
  height = 10,
  units = "cm",
  bg = "white",
  dpi = 300
)

# Save Figure 2 — effect sizes only
ggsave(
  filename = "output/figure_bacteria_richness_evenness_effect_only.tif",
  plot = figure_bacteria_richness_evenness_effect_only,
  width = 13.5,
  height = 10,
  units = "cm",
  bg = "white"
)
ggsave(
  filename = "output/figure_bacteria_richness_evenness_effect_only.png",
  plot = figure_bacteria_richness_evenness_effect_only,
  width = 13.5,
  height = 10,
  units = "cm",
  bg = "white",
  dpi = 300
)

#### FIGURE 3: Fungi richness & evenness ####

## Fungi richness plot
fungi_richness_plot <- ggplot() +
  geom_point(
    data = raw_data %>%
      filter(Response %in% c("Soil fungi richness", "Litter fungi richness")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi richness", "Litter fungi richness"))),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.5),
    alpha = 0.5
  ) +
  geom_errorbar(
    data = emmean %>%
      filter(Response %in% c("Soil fungi richness", "Litter fungi richness")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi richness", "Litter fungi richness"))),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>%
      filter(Response %in% c("Soil fungi richness", "Litter fungi richness"))%>%
      mutate(Response = factor(Response, levels = c("Soil fungi richness", "Litter fungi richness"))),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    plot.tag.position = c(0.03, 0.92)
  ) +
  scale_y_continuous(
    limits = c(
      0,
      raw_data %>%
        filter(Response %in% c("Soil fungi richness", "Litter fungi richness")) %>%
        pull(observed_value) %>%
        max()
    ),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = fungi_colours) +
  scale_fill_manual(values = fungi_colours) +
  labs(
    x = NULL,
    y = "Richness",
    tag = "(**a**)"
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil fungi richness" = "Soil",
      "Litter fungi richness" = "Litter"
    ))
  )

## Fungi richness effect plot
fungi_richness_density_nudged <- boot_densities %>%
  filter(Response %in% c("Soil fungi richness", "Litter fungi richness")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil fungi richness", "Litter fungi richness")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

fungi_richness_effect_nudged <- boot_effects %>%
  filter(Response %in% c("Soil fungi richness", "Litter fungi richness")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil fungi richness", "Litter fungi richness")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

fungi_richness_effect_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggridges::geom_density_ridges(
    data = fungi_richness_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type,
        group = interaction(regime, litter_type)),
    stat = "identity",
    scale = 0.75,
    colour = alpha("grey30", 0.5),
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = fungi_richness_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = fungi_richness_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  geom_richtext(
    data = intercept_annotations %>%
      filter(Response %in% c("Soil fungi richness", "Litter fungi richness")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi richness", "Litter fungi richness"))),
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    size = 3,
    fill = NA,
    label.color = NA
  ) +
  scale_y_continuous(
    breaks = unique(fungi_richness_density_nudged$regime_numeric),
    labels = unique(fungi_richness_density_nudged$regime),
    limits = c(0.5, max(fungi_richness_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(values = fungi_colours) +
  coord_flip() +
  my_theme +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(
    x = "Effect size",
    y = NULL
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil fungi richness" = "Soil",
      "Litter fungi richness" = "Litter"
    ))
  )

## Fungi evenness plot
fungi_evenness_plot <- ggplot() +
  geom_point(
    data = raw_data %>%
      filter(Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi evenness", "Litter fungi evenness"))),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.5),
    alpha = 0.5
  ) +
  geom_errorbar(
    data = emmean %>%
      filter(Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi evenness", "Litter fungi evenness"))),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>%
      filter(Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi evenness", "Litter fungi evenness"))),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_blank()
  ) +
  scale_y_continuous(
    limits = c(0.35, 1),
    breaks = scales::pretty_breaks(n = 4)
  ) +
  scale_color_manual(values = fungi_colours) +
  scale_fill_manual(values = fungi_colours) +
  labs(
    x = NULL,
    y = "Evenness",
    tag = "(**b**)"
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil fungi evenness" = "Soil",
      "Litter fungi evenness" = "Litter"
    ))
  )

## Fungi evenness effect plot
fungi_evenness_density_nudged <- boot_densities %>%
  filter(Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil fungi evenness", "Litter fungi evenness")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

fungi_evenness_effect_nudged <- boot_effects %>%
  filter(Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
  mutate(
    Response = factor(Response, levels = c("Soil fungi evenness", "Litter fungi evenness")),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

fungi_evenness_effect_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggridges::geom_density_ridges(
    data = fungi_evenness_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type,
        group = interaction(regime, litter_type)),
    stat = "identity",
    scale = 0.85,
    colour = alpha("grey30", 0.5),
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = fungi_evenness_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  geom_point(
    data = fungi_evenness_effect_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  geom_richtext(
    data = intercept_annotations %>%
      filter(Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
      mutate(Response = factor(Response, levels = c("Soil fungi evenness", "Litter fungi evenness"))),
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    size = 3,
    fill = NA,
    label.color = NA
  ) +
  scale_y_continuous(
    breaks = unique(fungi_evenness_density_nudged$regime_numeric),
    labels = unique(fungi_evenness_density_nudged$regime),
    limits = c(0.5, max(fungi_evenness_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(values = fungi_colours) +
  coord_flip() +
  my_theme +
  theme(strip.text = element_blank()) +
  labs(
    x = "Effect size",
    y = NULL
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil fungi evenness" = "Soil",
      "Litter fungi evenness" = "Litter"
    ))
  )

## Assemble Figure 3
figure_fungi_richness_evenness <- patchwork::wrap_plots(
  fungi_richness_plot, fungi_richness_effect_plot,
  fungi_evenness_plot, fungi_evenness_effect_plot,
  ncol = 1
)

figure_fungi_richness_evenness_final <- cowplot::plot_grid(
  figure_fungi_richness_evenness, dummy_legend_fungi, rel_widths = c(1, 0.24)
)

figure_fungi_richness_evenness_raw_only <- cowplot::plot_grid(
  patchwork::wrap_plots(
    fungi_richness_plot,
    fungi_evenness_plot +
      theme(axis.text.x = element_text(size = text_size)),
    ncol = 1
  ),
  dummy_legend_fungi,
  rel_widths = c(1, 0.24)
)

figure_fungi_richness_evenness_effect_only <- cowplot::plot_grid(
  patchwork::wrap_plots(
    fungi_richness_effect_plot +
      labs(tag = "(**a**)") +
      theme(
        strip.text = element_text(size = strip_size, face = "bold", vjust = 1),
        plot.tag.position = c(0.03, 0.97)
      ),
    fungi_evenness_effect_plot +
      labs(tag = "(**b**)"),
    ncol = 1
  ),
  dummy_legend_fungi,
  rel_widths = c(1, 0.24)
)

# Save Figure 3 — combined
ggsave(
  filename = "output/figure_fungi_richness_evenness.png",
  plot = figure_fungi_richness_evenness_final,
  width = 13.5,
  height = 20,
  units = "cm",
  bg = "white",
  dpi = 300
)
ggsave(
  filename = "output/figure_fungi_richness_evenness.tif",
  plot = figure_fungi_richness_evenness_final,
  width = 13.5,
  height = 20,
  units = "cm",
  bg = "white"
)

# Save Figure 3 — raw values and marginal means only
ggsave(
  filename = "output/figure_fungi_richness_evenness_raw_only.png",
  plot = figure_fungi_richness_evenness_raw_only,
  width = 13.5,
  height = 10,
  units = "cm",
  bg = "white",
  dpi = 300
)
ggsave(
  filename = "output/figure_fungi_richness_evenness_raw_only.tif",
  plot = figure_fungi_richness_evenness_raw_only,
  width = 13.5,
  height = 10,
  units = "cm",
  bg = "white"
)

# Save Figure 3 — effect sizes only
ggsave(
  filename = "output/figure_fungi_richness_evenness_effect_only.png",
  plot = figure_fungi_richness_evenness_effect_only,
  width = 13.5,
  height = 10,
  units = "cm",
  bg = "white",
  dpi = 300
)
ggsave(
  filename = "output/figure_fungi_richness_evenness_effect_only.tif",
  plot = figure_fungi_richness_evenness_effect_only,
  width = 13.5,
  height = 10,
  units = "cm",
  bg = "white"
)
