# Read in necessary libraries
require(ggtext)
require(tidyverse)

# Read in the data
load("data/alpha_diversity_plot_data.RData")

#### Shannon diversity plots wide ####

# Richness raw data plot
bacteria_diversity_plot <- ggplot() +
  geom_point(
    data = raw_data %>% filter(
      Response %in% c(
        "Soil bacteria shannon", "Litter bacteria shannon"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria shannon", "Litter bacteria shannon"
        ))
      ),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(
      jitter.width = 0.3,
      jitter.height = 0,
      dodge.width = 0.5
    ),
    alpha = 0.5,
  ) +
  geom_errorbar(
    data = emmean %>% filter(
      Response %in% c(
        "Soil bacteria shannon", "Litter bacteria shannon"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria shannon", "Litter bacteria shannon"
        ))
      ),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>% filter(
      Response %in% c(
        "Soil bacteria shannon", "Litter bacteria shannon"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria shannon", "Litter bacteria shannon"
        ))
      ),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    # Adjust the tag position to account for the strip text
    plot.tag.position = c(0.03, 0.92)
  ) +
  scale_y_continuous(
    limits = c(
      0, 
      raw_data %>%
        filter(Response %in% c(
          "Soil bacteria shannon", "Litter bacteria shannon",
          "Soil fungi shannon", "Litter fungi shannon"
        )) %>%
        pull(observed_value) %>%
        max()
    ),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_color_manual(
    values = litter_colours[1:2]
  ) +
  scale_fill_manual(
    values = litter_colours[1:2]
  ) +
  labs(
    x = NULL,
    y = "Shannon diversity",
    title = "Bacteria",
    tag = "(**a**)"
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil bacteria shannon" = "Soil",
      "Litter bacteria shannon" = "Litter"
    )))

# Display the plot
bacteria_diversity_plot

# Evenness raw data plot
fungi_diversity_plot <- ggplot() +
  geom_point(
    data = raw_data %>% filter(
      Response %in% c(
        "Soil fungi shannon", "Litter fungi shannon"
      )) %>%
      # Level the responses
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi shannon", "Litter fungi shannon"
        ))
      ),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(
      jitter.width = 0.3,
      jitter.height = 0,
      dodge.width = 0.5
    ),
    alpha = 0.5,
  ) +
  geom_errorbar(
    data = emmean %>% filter(
      Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
      # Add factor levels here too!
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi shannon", "Litter fungi shannon"
        ))
      ),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>% filter(
      Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
      # Add factor levels here too!
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi shannon", "Litter fungi shannon"
        ))
      ),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  ) +
  scale_y_continuous(
    limits = c(
      0, 
      raw_data %>%
        filter(Response %in% c(
          "Soil bacteria shannon", "Litter bacteria shannon",
          "Soil fungi shannon", "Litter fungi shannon"
        )) %>%
        pull(observed_value) %>%
        max()
    ),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_color_manual(
    values = litter_colours[1:2]
  ) +
  scale_fill_manual(
    values = litter_colours[1:2]
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Fungi",
  ) +
  facet_wrap(
    ~Response, ncol = 2,
    labeller = labeller(Response = c(
      "Soil fungi shannon" = "Soil",
      "Litter fungi shannon" = "Litter"
    )))

# Display the plot
fungi_diversity_plot

#### Shannon diversity plots long ####

# Richness raw data plot
bacteria_diversity_plot <- ggplot() +
  geom_point(
    data = raw_data %>% filter(
      Response %in% c(
        "Soil bacteria shannon", "Litter bacteria shannon"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria shannon", "Litter bacteria shannon"
        ))
      ),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(
      jitter.width = 0.3,
      jitter.height = 0,
      dodge.width = 0.5
    ),
    alpha = 0.5,
  ) +
  geom_errorbar(
    data = emmean %>% filter(
      Response %in% c(
        "Soil bacteria shannon", "Litter bacteria shannon"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria shannon", "Litter bacteria shannon"
        ))
      ),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>% filter(
      Response %in% c(
        "Soil bacteria shannon", "Litter bacteria shannon"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria shannon", "Litter bacteria shannon"
        ))
      ),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    # Adjust the tag position to account for the strip text
    plot.tag.position = c(0.03, 0.92)
  ) +
  scale_y_continuous(
    limits = c(
      0, 
      raw_data %>%
        filter(Response %in% c("Soil bacteria shannon", "Litter bacteria shannon")) %>%
        pull(observed_value) %>%
        max()
    ),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_color_manual(
    values = litter_colours[1:2]
  ) +
  scale_fill_manual(
    values = litter_colours[1:2]
  ) +
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
    )))

# Display the plot
bacteria_diversity_plot

# Evenness raw data plot
fungi_diversity_plot <- ggplot() +
  geom_point(
    data = raw_data %>% filter(
      Response %in% c(
        "Soil fungi shannon", "Litter fungi shannon"
      )) %>%
      # Level the responses
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi shannon", "Litter fungi shannon"
        ))
      ),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(
      jitter.width = 0.3,
      jitter.height = 0,
      dodge.width = 0.5
    ),
    alpha = 0.5,
  ) +
  geom_errorbar(
    data = emmean %>% filter(
      Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
      # Add factor levels here too!
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi shannon", "Litter fungi shannon"
        ))
      ),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>% filter(
      Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
      # Add factor levels here too!
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi shannon", "Litter fungi shannon"
        ))
      ),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(strip.text = element_blank()) +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(
    limits = c(
      0, 
      raw_data %>%
        filter(Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
        pull(observed_value) %>%
        max()
    ),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_color_manual(
    values = litter_colours[1:2]
  ) +
  scale_fill_manual(
    values = litter_colours[1:2]
  ) +
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
    )))

# Display the plot
fungi_diversity_plot

#### bacteria effect richness & evenness plots ####

# Density data for richness
richness_effect_density_nudged <- effect_density %>% 
  filter(Response %in% c("Soil bacteria shannon", "Litter bacteria shannon")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil bacteria shannon", "Litter bacteria shannon"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Density data for evenness
evenness_effect_density_nudged <- effect_density %>%
  filter(Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil fungi shannon", "Litter fungi shannon"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Mean data for richness
richness_effect_mean_nudged <- effect_mean %>%
  filter(Response %in% c("Soil bacteria shannon", "Litter bacteria shannon")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil bacteria shannon", "Litter bacteria shannon"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Mean data for evenness
evenness_effect_mean_nudged <- effect_mean %>%
  filter(Response %in% c("Soil fungi shannon", "Litter fungi shannon")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil fungi shannon", "Litter fungi shannon"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Effect plot for richness
bacteria_diversity_effect_plot <- ggplot() +
  # Baseline (the grand mean)
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  ggridges::geom_density_ridges(
    data = richness_effect_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type, 
        group = interaction(regime, litter_type)),
    stat = "identity", 
    scale = 0.75, 
    colour = alpha("grey30", 0.5),
    alpha = 0.7
  ) +
  # Add error bars
  geom_errorbar(
    data = richness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, 
        xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  # Add points
  geom_point(
    data = richness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = unique(richness_effect_density_nudged$regime_numeric),
    labels = unique(richness_effect_density_nudged$regime),
    limits = c(0.5, max(richness_effect_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(
    values = litter_colours[1:2],
  ) +
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
    )))

# Effect plot for evenness
fungi_diversity_effect_plot <- ggplot() +
  # Baseline (the grand mean)
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  ggridges::geom_density_ridges(
    data = evenness_effect_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type, 
        group = interaction(regime, litter_type)),
    stat = "identity", 
    scale = 0.75, 
    colour = alpha("grey30", 0.5),
    alpha = 0.7
  ) +
  # Add error bars
  geom_errorbar(
    data = evenness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, 
        xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  # Add points
  geom_point(
    data = evenness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = unique(evenness_effect_density_nudged$regime_numeric),
    labels = unique(evenness_effect_density_nudged$regime),
    limits = c(0.5, max(evenness_effect_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(
    values = litter_colours[1:2],
  ) +
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
    )))

figure_bacteria_diversity_evenness <-patchwork::wrap_plots(
  bacteria_diversity_plot, bacteria_diversity_effect_plot,
  fungi_diversity_plot, fungi_diversity_effect_plot,
  ncol = 1
)

figure_bacteria_diversity_evenness_final <- cowplot::plot_grid(
  figure_bacteria_diversity_evenness, dummy_legend, rel_widths = c(1, 0.24)
)

# Save the final figure
ggsave(
  filename = "output/figure_shannon_diversity_long.pdf",
  plot = figure_bacteria_diversity_evenness_final,
  width = 14,
  height = 20, 
  units = "cm"
)

#### Bacteria mean richness & evenness plots ####

# Richness raw data plot
bacteria_richness_plot <- ggplot() +
  geom_point(
    data = raw_data %>% filter(
      Response %in% c(
        "Soil bacteria richness", "Litter bacteria richness"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria richness", "Litter bacteria richness"
        ))
      ),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(
      jitter.width = 0.3,
      jitter.height = 0,
      dodge.width = 0.5
    ),
    alpha = 0.5,
  ) +
  geom_errorbar(
    data = emmean %>% filter(
      Response %in% c(
        "Soil bacteria richness", "Litter bacteria richness"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria richness", "Litter bacteria richness"
        ))
      ),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>% filter(
      Response %in% c(
        "Soil bacteria richness", "Litter bacteria richness"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria richness", "Litter bacteria richness"
        ))
      ),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    # Adjust the tag position to account for the strip text
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
  scale_color_manual(
    values = litter_colours[1:2]
  ) +
  scale_fill_manual(
    values = litter_colours[1:2]
  ) +
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
    )))

# Display the plot
bacteria_richness_plot

# Evenness raw data plot
bacteria_evenness_plot <- ggplot() +
  geom_point(
    data = raw_data %>% filter(
      Response %in% c(
        "Soil bacteria evenness", "Litter bacteria evenness"
      )) %>%
      # Level the responses
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria evenness", "Litter bacteria evenness"
        ))
      ),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(
      jitter.width = 0.3,
      jitter.height = 0,
      dodge.width = 0.5
    ),
    alpha = 0.5,
  ) +
  geom_errorbar(
    data = emmean %>% filter(
      Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
      # Add factor levels here too!
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria evenness", "Litter bacteria evenness"
        ))
      ),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>% filter(
      Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
      # Add factor levels here too!
      mutate(
        Response = factor(Response, levels = c(
          "Soil bacteria evenness", "Litter bacteria evenness"
        ))
      ),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(strip.text = element_blank()) +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(
    limits = c(
      0, 
      raw_data %>%
        filter(Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
        pull(observed_value) %>%
        max()
    ),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_color_manual(
    values = litter_colours[1:2]
  ) +
  scale_fill_manual(
    values = litter_colours[1:2]
  ) +
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
    )))

# Display the plot
bacteria_evenness_plot

#### Bacteria effect richness & evenness plots ####

# Density data for richness
richness_effect_density_nudged <- effect_density %>% 
  filter(Response %in% c("Soil bacteria richness", "Litter bacteria richness")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil bacteria richness", "Litter bacteria richness"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Density data for evenness
evenness_effect_density_nudged <- effect_density %>%
  filter(Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil bacteria evenness", "Litter bacteria evenness"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Mean data for richness
richness_effect_mean_nudged <- effect_mean %>%
  filter(Response %in% c("Soil bacteria richness", "Litter bacteria richness")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil bacteria richness", "Litter bacteria richness"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Mean data for evenness
evenness_effect_mean_nudged <- effect_mean %>%
  filter(Response %in% c("Soil bacteria evenness", "Litter bacteria evenness")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil bacteria evenness", "Litter bacteria evenness"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Effect plot for richness
bacteria_richness_effect_plot <- ggplot() +
  # Baseline (the grand mean)
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  ggridges::geom_density_ridges(
    data = richness_effect_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type, 
        group = interaction(regime, litter_type)),
    stat = "identity", 
    scale = 0.75, 
    colour = alpha("grey30", 0.5),
    alpha = 0.7
  ) +
  # Add error bars
  geom_errorbar(
    data = richness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, 
        xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  # Add points
  geom_point(
    data = richness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = unique(richness_effect_density_nudged$regime_numeric),
    labels = unique(richness_effect_density_nudged$regime),
    limits = c(0.5, max(richness_effect_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(
    values = litter_colours[1:2],
  ) +
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
    )))

# Effect plot for evenness
bacteria_evenness_effect_plot <- ggplot() +
  # Baseline (the grand mean)
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  ggridges::geom_density_ridges(
    data = evenness_effect_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type, 
        group = interaction(regime, litter_type)),
    stat = "identity", 
    scale = 0.85, 
    colour = alpha("grey30", 0.5),
    alpha = 0.7
  ) +
  # Add error bars
  geom_errorbar(
    data = evenness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, 
        xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  # Add points
  geom_point(
    data = evenness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = unique(evenness_effect_density_nudged$regime_numeric),
    labels = unique(evenness_effect_density_nudged$regime),
    limits = c(0.5, max(evenness_effect_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(
    values = litter_colours[1:2],
  ) +
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
    )))

figure_bacteria_richness_evenness <-patchwork::wrap_plots(
  bacteria_richness_plot, bacteria_richness_effect_plot,
  bacteria_evenness_plot, bacteria_evenness_effect_plot,
  ncol = 1
)

figure_bacteria_richness_evenness_final <- cowplot::plot_grid(
  figure_bacteria_richness_evenness, dummy_legend, rel_widths = c(1, 0.24)
)

# Save the final figure
ggsave(
  filename = "output/figure_bacteria_richness_evenness.pdf",
  plot = figure_bacteria_richness_evenness_final,
  width = 14,
  height = 20, 
  units = "cm"
)

#### Fungi mean richness & evenness plots ####

# Richness raw data plot
fungi_richness_plot <- ggplot() +
  geom_point(
    data = raw_data %>% filter(
      Response %in% c(
        "Soil fungi richness", "Litter fungi richness"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi richness", "Litter fungi richness"
        ))
      ),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(
      jitter.width = 0.3,
      jitter.height = 0,
      dodge.width = 0.5
    ),
    alpha = 0.5,
  ) +
  geom_errorbar(
    data = emmean %>% filter(
      Response %in% c(
        "Soil fungi richness", "Litter fungi richness"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi richness", "Litter fungi richness"
        ))
      ),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>% filter(
      Response %in% c(
        "Soil fungi richness", "Litter fungi richness"
      )) %>%
      # Level the resonses
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi richness", "Litter fungi richness"
        ))
      ),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    # Adjust the tag position to account for the strip text
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
  scale_color_manual(
    values = litter_colours[1:2]
  ) +
  scale_fill_manual(
    values = litter_colours[1:2]
  ) +
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
    )))

# Display the plot
fungi_richness_plot

# Evenness raw data plot
fungi_evenness_plot <- ggplot() +
  geom_point(
    data = raw_data %>% filter(
      Response %in% c(
        "Soil fungi evenness", "Litter fungi evenness"
      )) %>%
      # Level the responses
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi evenness", "Litter fungi evenness"
        ))
      ),
    aes(x = regime, y = observed_value, color = litter_type),
    position = position_jitterdodge(
      jitter.width = 0.3,
      jitter.height = 0,
      dodge.width = 0.5
    ),
    alpha = 0.5,
  ) +
  geom_errorbar(
    data = emmean %>% filter(
      Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
      # Add factor levels here too!
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi evenness", "Litter fungi evenness"
        ))
      ),
    aes(x = regime, y = emmean, ymin = lower_CI, ymax = upper_CI, group = litter_type),
    colour = "black",
    position = position_dodge(width = 0.5),
    width = 0.1,
    size = 1
  ) +
  geom_point(
    data = emmean %>% filter(
      Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
      # Add factor levels here too!
      mutate(
        Response = factor(Response, levels = c(
          "Soil fungi evenness", "Litter fungi evenness"
        ))
      ),
    aes(x = regime, y = emmean, fill = litter_type),
    shape = 21,
    colour = "black",
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  my_theme +
  theme(strip.text = element_blank()) +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(
    limits = c(
      0, 
      raw_data %>%
        filter(Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
        pull(observed_value) %>%
        max()
    ),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_color_manual(
    values = litter_colours[1:2]
  ) +
  scale_fill_manual(
    values = litter_colours[1:2]
  ) +
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
    )))

# Display the plot
fungi_evenness_plot

#### Fungi effect richness & evenness plots ####

# Density data for richness
richness_effect_density_nudged <- effect_density %>% 
  filter(Response %in% c("Soil fungi richness", "Litter fungi richness")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil fungi richness", "Litter fungi richness"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Density data for evenness
evenness_effect_density_nudged <- effect_density %>%
  filter(Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil fungi evenness", "Litter fungi evenness"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Mean data for richness
richness_effect_mean_nudged <- effect_mean %>%
  filter(Response %in% c("Soil fungi richness", "Litter fungi richness")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil fungi richness", "Litter fungi richness"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Mean data for evenness
evenness_effect_mean_nudged <- effect_mean %>%
  filter(Response %in% c("Soil fungi evenness", "Litter fungi evenness")) %>%
  mutate(
    Response = factor(Response, levels = c(
      "Soil fungi evenness", "Litter fungi evenness"
    )),
    regime_numeric = as.numeric(factor(regime)),
    regime_nudged = case_when(
      litter_type == "Broadleaf" ~ regime_numeric - 0.15,
      litter_type == "Grass" ~ regime_numeric + 0.15,
      TRUE ~ regime_numeric
    )
  )

# Effect plot for richness
fungi_richness_effect_plot <- ggplot() +
  # Baseline (the grand mean)
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  ggridges::geom_density_ridges(
    data = richness_effect_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type, 
        group = interaction(regime, litter_type)),
    stat = "identity", 
    scale = 0.75, 
    colour = alpha("grey30", 0.5),
    alpha = 0.7
  ) +
  # Add error bars
  geom_errorbar(
    data = richness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, 
        xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  # Add points
  geom_point(
    data = richness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = unique(richness_effect_density_nudged$regime_numeric),
    labels = unique(richness_effect_density_nudged$regime),
    limits = c(0.5, max(richness_effect_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(
    values = litter_colours[1:2],
  ) +
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
    )))

# Effect plot for evenness
fungi_evenness_effect_plot <- ggplot() +
  # Baseline (the grand mean)
  geom_vline(
    xintercept = 0,
    linetype = "dotted"
  ) +
  ggridges::geom_density_ridges(
    data = evenness_effect_density_nudged,
    aes(x = Effect, y = regime_nudged, height = Density, fill = litter_type, 
        group = interaction(regime, litter_type)),
    stat = "identity", 
    scale = 0.85, 
    colour = alpha("grey30", 0.5),
    alpha = 0.7
  ) +
  # Add error bars
  geom_errorbar(
    data = evenness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, 
        xmin = lower_ci, xmax = upper_ci, group = litter_type),
    color = "black",
    width = 0,
    linewidth = 1
  ) +
  # Add points
  geom_point(
    data = evenness_effect_mean_nudged,
    aes(x = effect_mean, y = regime_nudged, fill = litter_type),
    shape = 21,
    size = 3,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = unique(evenness_effect_density_nudged$regime_numeric),
    labels = unique(evenness_effect_density_nudged$regime),
    limits = c(0.5, max(evenness_effect_density_nudged$regime_numeric) + 0.5)
  ) +
  scale_fill_manual(
    values = litter_colours[1:2],
  ) +
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
    )))

figure_fungi_richness_evenness <-patchwork::wrap_plots(
  fungi_richness_plot, fungi_richness_effect_plot,
  fungi_evenness_plot, fungi_evenness_effect_plot,
  ncol = 1
)

figure_fungi_richness_evenness_final <- cowplot::plot_grid(
  figure_fungi_richness_evenness, dummy_legend, rel_widths = c(1, 0.24)
)

# Save the final figure
ggsave(
  filename = "output/figure_fungi_richness_evenness.pdf",
  plot = figure_fungi_richness_evenness_final,
  width = 14,
  height = 20, 
  units = "cm"
)
