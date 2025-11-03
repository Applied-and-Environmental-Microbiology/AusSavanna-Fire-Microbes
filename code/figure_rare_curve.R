# Read in necessary libraries
require(ggtext)
require(patchwork)
require(tidyverse)

# Read in the data
load("data/rare_curve_bacteria.RData")
load("data/rare_curve_fungi.RData")

# Max depth for the shared x axis
max_depth <- max(
  max_depth_bacteria,
  max_depth_fungi
)

# Plot the results
rare_curves <- wrap_plots(
  
  # Bactaria
  map_dfr(rare_curve_bacteria, bind_rows) %>%
    bind_cols(sample_id = rownames(otus_trans_bacteria),.) %>%
    # Organise the data
    pivot_longer(-sample_id) %>%
    drop_na() %>%
    mutate(depth = as.numeric(str_replace(name, 'N', ''))) %>%
    select(-name) %>%
    # Plot with ggplot2
    ggplot(aes(x = depth, y = value, sample_id = sample_id)) +
    geom_vline(xintercept = min_depth_bacteria, colour = 'grey') +
    geom_line() +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
      axis.text.x = element_blank(),
      axis.text.y = element_text(colour = 'black'),
      aspect.ratio = 1,
      plot.tag = element_markdown(),
    ) + 
    scale_x_continuous(
      labels = scales::comma_format(),
      limits = c(0, max_depth)
    ) +
    scale_y_continuous(
      labels = scales::comma_format()
    ) +
    labs(
      tag = '(**a**)',
      y = 'Number of ASVs',
      x = NULL
    ),
  
  # Fungi
  map_dfr(rare_curve_fungi, bind_rows) %>%
    bind_cols(sample_id = rownames(otus_trans_fungi),.) %>%
    # Organise the data
    pivot_longer(-sample_id) %>%
    drop_na() %>%
    mutate(depth = as.numeric(str_replace(name, 'N', ''))) %>%
    select(-name) %>%
    # Plot with ggplot2
    ggplot(aes(x = depth, y = value, sample_id = sample_id)) +
    geom_vline(xintercept = min_depth_fungi, colour = 'grey') +
    geom_line() +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
      axis.text.x = element_text(colour = 'black'),
      axis.text.y = element_text(colour = 'black'),
      aspect.ratio = 1,
      plot.tag = element_markdown(),
    ) + 
    scale_x_continuous(
      labels = scales::comma_format(),
      limits = c(0, max_depth)
    ) +
    scale_y_continuous(
      labels = scales::comma_format()
    ) +
    labs(
      tag = '(**b**)',
      y = 'Number of OTUs',
      x = 'Number of reads'
    ),
  
  # Number of columns
  ncol = 1

)

# View the plot
rare_curves

# Save the plot
ggsave(
  filename = 'output/rare_curves.png',
  plot = rare_curves,
  width = 10,
  height = 20,
  units = 'cm',
  dpi = 300
)
