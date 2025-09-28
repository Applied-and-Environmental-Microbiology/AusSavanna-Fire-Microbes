
# Required packages and functions
require(Biostrings)
require(vegan)
require(SRS)
require(data.table)
require(tidyverse)

# Helper function for computing shannon diversity
shannon_div <- function(x) {
  p <- x / sum(x)
  -sum(p * log(p))
}

# Helper function for mapping from sequences to OTU IDs
sequence_to_otu_mapping <- data.frame(
  sequence = as.character(sequences),
  otu_id = names(sequences),
  stringsAsFactors = FALSE
)

# Read in the raw data files: OTU IDs in the OTU and taxonomy tables need to be
# converted from actual sequence to OTU IDs (names in the fatsa file)
otus_raw <- fread("data/ASVs_abundance_filtered.txt")
taxonomy_raw <- fread("data/taxonomy_abundance_filtered.txt")
sequences <- readDNAStringSet("data/ASVs.fasta")
metadata <- fread("data/metadata.csv") %>%
  select(
    sample_id = SampleID, regime = Regime, block = Block, plot = Plot,
    quadrat = Quadrat, litter_type = Veg, community = Community
  )

# Replace sequence IDs with OTU IDs in both tables
otus <- otus_raw %>%
  left_join(sequence_to_otu_mapping, by = c("ASV_ID" = "sequence")) %>%
  select(ASV_ID = otu_id, everything(), -ASV_ID) %>%  # Keep ASV_ID column name for compatibility
  filter(!is.na(ASV_ID))  # Remove any sequences not found in the mapping

taxonomy <- taxonomy_raw %>%
  left_join(sequence_to_otu_mapping, by = c("ASV_ID" = "sequence")) %>%
  select(ASV_ID = otu_id, everything(), -ASV_ID) %>%  # Keep ASV_ID column name for compatibility  
  filter(!is.na(ASV_ID))  


# (1) Assess minimum read depth ------------------------------------------------

# Create a long format tibble of OTUs
otus_tibble <- otus %>%
  # Transpose the OTU table
  column_to_rownames(var = 'ASV_ID') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_id') %>%
  as_tibble() %>%
  # Convert to long format
  pivot_longer(-sample_id) %>%
  # Remove zero abundance OTUs
  filter(value != 0)

# Compute the read depth per sample_id
sample_id_depth <- otus_tibble %>%
  # Count total reads per sample_id
  group_by(sample_id) %>%
  summarise(n_seqs = sum(value)) %>%
  arrange(n_seqs) %>%
  print(n = 50)

# 4,300 appears to be a reasonable lowest read depth

# Evaluate the range of read depth
message("Range of read depth: ", min(sample_id_depth$n_seqs), " - ", max(sample_id_depth$n_seqs))

# Visualise sample_id depth and range
sample_id_depth %>%
  ggplot(aes(x = 1:nrow(.), y = n_seqs)) +
  geom_line() +
  geom_point()
sample_id_depth %>%
  ggplot(aes(x = 1, y = n_seqs)) +
  geom_jitter()

# Set the minimum read depth threshold
min_depth <- 4300
max_depth <- max(sample_id_depth$n_seqs)

# Create a vector of low abundance sample_ids
low_abundance_sample_ids <- sample_id_depth %>%
  filter(n_seqs < min_depth) %>%
  pull(sample_id)

# OTU richness per sample_id
otu_richness = otus_tibble %>%
  group_by(sample_id) %>%
  summarise(n_otus = n_distinct(name)) %>%
  arrange(n_otus)

# Range of OTU richness
message("Range of OTU richness: ", min(otu_richness$n_otus), " - ", max(otu_richness$n_otus))

# Evaluate the relationship between read depth and OTU richness
sample_id_depth %>%
  left_join(otu_richness, by = "sample_id") %>%
  # Remove the low abundance sample_ids
  filter(n_seqs >= min_depth) %>%
  ggplot(aes(x = n_seqs, y = n_otus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(
    # Add the r2 label
    aes(label = after_stat(rr.label))
  ) +
  labs(x = "Read depth", y = "OTU richness")

# The relationship is weak (r2 = 0.03), suggesting that the sequencing depth
# does not strongly influences the OTU richness

# (2) Species accumulation curves ----------------------------------------------

# Vegan expects sample_id in rows and OTU in columns
otus_trans <- otus %>%
  # Remove the low abundance sample_ids
  select(-all_of(low_abundance_sample_ids)) %>%
  # Transpose the OTU table
  column_to_rownames(var = 'ASV_ID') %>%
  t() %>%
  as.data.frame()

# Create a species accumulation curve
sp_curve <- specaccum(otus_trans, 'random')
plot(sp_curve)

# Check out another approach with multiple options
sp_curve_2 <- poolaccum(otus_trans)
plot(sp_curve_2)

# Organise the data for plotting
sp_curve_data <- data.frame(
  summary(sp_curve_2)$S, check.names = F
) %>%
  # Select relevant columns
  select(
    N = N, S = S, 
    lowCI = '2.5%', upCI = '97.5%', 
    SD = Std.Dev
  )

# Plot the species accumulation curve with ggplot2
plot_sp_curve <- ggplot(
  data = sp_curve_data, 
  aes(x = N, y = S)
) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = lowCI,
                  ymax = upCI),
              alpha = 0.5,
              colour = "gray70") +
  # Add observed richness line 
  geom_line() +
  scale_y_continuous(
    limits = c(0, max(sp_curve_data$upCI)),
    labels = scales::comma_format(big.mark = ','),
    breaks = scales::pretty_breaks(n = 5),
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    axis.ticks = element_blank(),
    aspect.ratio = 1
  ) + 
  labs(x = 'Number of sample_ids',
       y = 'Number of OTUs')

# Display the plot
print(plot_sp_curve)

# (3) Rarefaction curves -------------------------------------------------------

# Create a rarefaction curve
rare_curve = rarecurve(otus_trans, step = 100)

# Plot the rarefaction curves with ggplot2
plot_rare_curve <- map_dfr(rare_curve, bind_rows) %>%
  bind_cols(sample_id = rownames(otus_trans),.) %>%
  # Organise the data
  pivot_longer(-sample_id) %>%
  drop_na() %>%
  mutate(depth = as.numeric(str_replace(name, 'N', ''))) %>%
  select(-name) %>%
  # Plot with ggplot2
  ggplot(aes(x = depth, y = value, sample_id = sample_id)) +
  geom_vline(xintercept = min_depth, colour = 'grey') +
  geom_line() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    axis.ticks = element_blank(),
    aspect.ratio = 1
  ) + 
  scale_x_continuous(
    labels = scales::comma_format(),
    limits = c(0, max_depth)
  ) +
  scale_y_continuous(
    labels = scales::comma_format()
  ) +
  xlab('Number of reads') +
  ylab('Number of OTUs')

# Display the plot
print(plot_rare_curve)

# (4) Compute alpha-diversity metrics ------------------------------------------

# For alpha-diversity analyses we will rarefy to 4,300 reads (minimum depth)
# using the SRS method

# Remove the low abundance sample_ids from the OTU table
otus_2 <- otus %>%
  # Remove the low abundance sample_ids
  select(-all_of(low_abundance_sample_ids))

# SRS expects columns are sample_ids and rows are OTUs
otus_srs <- SRS(
  # OTU ID's as rownames
  data = otus_2 %>% column_to_rownames(var = 'ASV_ID'),
  Cmin = min_depth,  # <-- change to the minimum depth of the prevalence filtered OTU table for beta-diversity analyses
  set_seed = TRUE,  # <-- for reproducibility
  seed = 1986
) %>%
  # Convert the SRS output to a data frame
  as_tibble() %>%
  # Bind the OTU_ID's and SRS normalised counts
  bind_cols(
    tibble(otu_id = otus_2[["ASV_ID"]]),
    .
  ) %>%
  # Remove any OTUs that are now zero abundance
  filter(rowSums(select(., -otu_id)) > 0)

# Compute alpha-diversity metrics for rarefied OTUs
alpha_diversity_bacteria = otus_srs %>%
  # Convert to long format
  pivot_longer(
    -otu_id,
    names_to = "sample_id",
    values_to = "count"
  ) %>%
  # Remove zero abundance OTUs
  filter(count != 0) %>%
  group_by(sample_id) %>%
  summarise(
    richness_bacteria = n_distinct(otu_id),
    shannon_bacteria = shannon_div(count)
  ) %>%
  mutate(
    evenness_bacteria = shannon_bacteria / log(richness_bacteria)
  )

# Join the alpha-diversity metrics with the metadata
metadata_div <- metadata %>%
  inner_join(alpha_diversity_bacteria, by = "sample_id")

# (5) Save the outputs ---------------------------------------------------------

# Save the abundance filtered OTU table
fwrite(
  otus_srs,
  "output/otus.txt", 
  sep = "\t"
)

# Save the prevalence filtered (for use in gemelli) - these should not be normalised or rarefied
# I need to generate two files, one for "Leaf" and one for "Soil"
leaf_samples <- metadata %>%
  filter(community == "Leaf") %>%
  pull(sample_id)
soil_samples <- metadata %>%
  filter(community == "Soil") %>%
  pull(sample_id)

# Save the leaf prevalence filtered OTU table
fwrite(
  otus %>%
    # Select only leaf or soil sample_ids
    select(otu_id = ASV_ID, all_of(leaf_samples)) %>%
    # Remove the low abundance sample_ids
    select(-all_of(low_abundance_sample_ids)) %>%
    filter(rowSums(select(., -otu_id) > 0) >= 2),  # <-- at least 2 sample_ids
  "output/otus_leaf_prevalence_filtered.txt", 
  sep = "\t"
)
fwrite(
  otus %>%
    # Select only leaf or soil sample_ids
    select(otu_id = ASV_ID, all_of(soil_samples)) %>%
    filter(rowSums(select(., -otu_id) > 0) >= 2),  # <-- at least 2 sample_ids
  "output/otus_soil_prevalence_filtered.txt", 
  sep = "\t"
)

# Save the taxonomy table
fwrite(
  taxonomy %>%
    dplyr::rename(otu_id = ASV_ID) %>%
    # Remove otus not in the normalised table for consistency
    filter(otu_id %in% otus_srs$otu_id),
  "output/taxonomy.txt", 
  sep = "\t"
)

# Extract OTU IDs from sequence names (everything before the first space)
sequence_otu_ids <- sub(" .*", "", names(sequences))

# Filter sequences where the extracted OTU ID matches your list
filtered_sequences <- sequences[sequence_otu_ids %in% otus_srs$otu_id]

# Save the filtered sequences to a new FASTA file
writeXStringSet(
  filtered_sequences,
  filepath = "output/sequences.fasta"
)

# (6) Compute distance ---------------------------------------------------------

# Create a directory for gemelli output
system("mkdir -p data/distances_leaf")
system("mkdir -p data/distances_soil")

# Convert to biom format_leaf
system('conda run -n gemelli_env biom convert -i output/otus_leaf_prevalence_filtered.txt -o data/distances_leaf/otus.biom --table-type="OTU table" --to-json')
system('conda run -n gemelli_env biom convert -i output/otus_soil_prevalence_filtered.txt -o data/distances_soil/otus.biom --table-type="OTU table" --to-json')

# Run gemelli
system("conda run -n gemelli_env gemelli rpca --in-biom data/distances_leaf/otus.biom --output-dir data/distances_leaf/")
system("conda run -n gemelli_env gemelli rpca --in-biom data/distances_soil/otus.biom --output-dir data/distances_soil/")

# List the output files
system("ls data/distances_leaf")
system("ls data/distances_soil")

# Read the ordination file
ordination_lines_leaf <- readLines("data/distances_leaf/ordination.txt")
ordination_lines_soil <- readLines("data/distances_soil/ordination.txt")

# Find the line indices for different sections
proportion_line_leaf <- which(grepl("Proportion explained", ordination_lines_leaf))[1]  # Get first match
proportion_line_soil <- which(grepl("Proportion explained", ordination_lines_soil))[1]  # Get first match

site_line_leaf <- which(grepl("^Site", ordination_lines_leaf))[1]  # Get first match, use ^ for start of line
site_line_soil <- which(grepl("^Site", ordination_lines_soil))[1]  # Get first match, use ^ for start of line

# Extract the proportion explained values
proportions_leaf <- as.numeric(strsplit(ordination_lines_leaf[proportion_line_leaf + 1], "\t")[[1]])
proportions_soil <- as.numeric(strsplit(ordination_lines_soil[proportion_line_soil + 1], "\t")[[1]])

prop1_leaf <- round(proportions_leaf[1] * 100, 0)  # Convert to percentage and round
prop2_leaf <- round(proportions_leaf[2] * 100, 0)
prop1_soil <- round(proportions_soil[1] * 100, 0)  # Convert to percentage and round
prop2_soil <- round(proportions_soil[2] * 100, 0)

# Read the site coordinates (starting from the line after "Site")
site_start_leaf <- site_line_leaf + 1
site_end_leaf <- length(ordination_lines_leaf)
site_start_soil <- site_line_soil + 1
site_end_soil <- length(ordination_lines_soil)

# Find where site data ends (look for empty lines or "Biplot"/"Site constraints")
for (i in site_start_leaf:length(ordination_lines_leaf)) {
  if (ordination_lines_leaf[i] == "" || grepl("^Biplot", ordination_lines_leaf[i]) || grepl("^Site constraints", ordination_lines_leaf[i])) {
    site_end <- i - 1
    break
  }
}
for (i in site_start_soil:length(ordination_lines_soil)) {
  if (ordination_lines_soil[i] == "" || grepl("^Biplot", ordination_lines_soil[i]) || grepl("^Site constraints", ordination_lines_soil[i])) {
    site_end <- i - 1
    break
  }
}

site_data_leaf <- ordination_lines_leaf[site_start_leaf:site_end_leaf]
site_data_soil <- ordination_lines_soil[site_start_soil:site_end_soil]

# Parse site data into a data frame
site_coords_leaf <- data.frame(
  sample_id = character(),
  PC1 = numeric(),
  PC2 = numeric(),
  stringsAsFactors = FALSE
)
site_coords_soil <- data.frame(
  sample_id = character(),
  PC1 = numeric(),
  PC2 = numeric(),
  stringsAsFactors = FALSE
)

for (line in site_data_leaf) {
  if (line != "" && !grepl("^\t", line)) {  # Skip empty lines and lines starting with tab
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      site_coords_leaf <- rbind(site_coords_leaf, data.frame(
        sample_id = parts[1],
        PC1 = as.numeric(parts[2]),
        PC2 = as.numeric(parts[3]),
        stringsAsFactors = FALSE
      )) %>%
        filter(!sample_id %in% c("Biplot", "Site constraints"))
    }
  }
}
for (line in site_data_soil) {
  if (line != "" && !grepl("^\t", line)) {  # Skip empty lines and lines starting with tab
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      site_coords_soil <- rbind(site_coords_soil, data.frame(
        sample_id = parts[1],
        PC1 = as.numeric(parts[2]),
        PC2 = as.numeric(parts[3]),
        stringsAsFactors = FALSE
      )) %>%
        filter(!sample_id %in% c("Biplot", "Site constraints"))
    }
  }
} 

# Create the final data frame with your specified column names
ordination_results_leaf <- site_coords_leaf %>%
  dplyr::rename(
    !!paste0("pco1_", prop1_leaf, "_leaf_bacteria") := PC1,
    !!paste0("pco2_", prop2_leaf, "_leaf_bacteria") := PC2
  ) %>%
  select(sample_id, everything())
ordination_results_soil <- site_coords_soil %>%
  dplyr::rename(
    !!paste0("pco1_", prop1_soil, "_soil_bacteria") := PC1,
    !!paste0("pco2_", prop2_soil, "_soil_bacteria") := PC2
  ) %>%
  select(sample_id, everything())

# Add the ordination results to the metadata
metadata_beta_div <- metadata_div %>%
  full_join(ordination_results_leaf, by = "sample_id") %>%
  full_join(ordination_results_soil, by = "sample_id")

# (7) Save or move results for statistical analyses ----------------------------

# Save the data for statistical analyses
fwrite(
  metadata_beta_div,
  "../../data/sample_metadata_bacteria.txt", 
  sep = "\t"
)

# Move the distance matrices to the main data directory
system("mv data/distances_leaf/distance-matrix.tsv ../../data/distance_bacteria_leaf.txt")
system("mv data/distances_soil/distance-matrix.tsv ../../data/distance_bacteria_soil.txt")
