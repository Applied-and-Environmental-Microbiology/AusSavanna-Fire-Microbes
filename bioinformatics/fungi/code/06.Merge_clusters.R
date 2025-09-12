
# Set worke from the R project root

# Load packages
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(tidyverse))

# Clusters across ranks, rename pseudo taxa, select reference IDs ##############
taxonomy_1 <- fread("code/tmp_clusters/species_clusters.txt") %>%
  select(OTU_ID, reference_ID, species, cutoff, confidence, score) %>%
  inner_join(
    .,
    fread("code/tmp_clusters/genus_clusters.txt") %>%
      select(OTU_ID, genus),
    by = "OTU_ID"
  ) %>%
  inner_join(
    .,
    fread("code/tmp_clusters/family_clusters.txt") %>%
      select(OTU_ID, family),
    by = "OTU_ID"
  ) %>%
  inner_join(
    .,
    fread("code/tmp_clusters/order_clusters.txt") %>%
      select(OTU_ID, order),
    by = "OTU_ID"
  ) %>%
  inner_join(
    .,
    fread("code/tmp_clusters/class_clusters.txt") %>%
      select(OTU_ID, class),
    by = "OTU_ID"
  ) %>%
  inner_join(
    .,
    fread("code/tmp_clusters/phylum_clusters.txt") %>%
      select(OTU_ID, phylum),
    by = "OTU_ID"
  ) %>%
  select(
    OTU_ID, reference_ID, phylum, class, order, family, genus, species, score, cutoff
  ) %>%
  # Rename pseudo taxa
  group_by(phylum) %>%
  mutate(
    phylum = case_when(
      str_detect(phylum, "_pseudo_") ~ paste0(str_extract(phylum, "^[^_]+"), "_pseudo_phylum_", sprintf("%04d", cur_group_id())),
      TRUE ~ phylum
    )) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(
    class = case_when(
      str_detect(class, "_pseudo_") ~ paste0(str_extract(class, "^[^_]+"), "_pseudo_class_", sprintf("%04d", cur_group_id())),
      TRUE ~ class
    )) %>%
  ungroup() %>%
  group_by(order) %>%
  mutate(
    order = case_when(
      str_detect(order, "_pseudo_") ~ paste0(str_extract(order, "^[^_]+"), "_pseudo_order_", sprintf("%04d", cur_group_id())),
      TRUE ~ order
    )) %>%
  ungroup() %>%
  group_by(family) %>%
  mutate(
    family = case_when(
      str_detect(family, "_pseudo_") ~ paste0(str_extract(family, "^[^_]+"), "_pseudo_family_", sprintf("%04d", cur_group_id())),
      TRUE ~ family
    )) %>%
  ungroup() %>%
  group_by(genus) %>%
  mutate(
    genus = case_when(
      str_detect(genus, "_pseudo_") ~ paste0(str_extract(genus, "^[^_]+"), "_pseudo_genus_", sprintf("%04d", cur_group_id())),
      TRUE ~ genus
    )) %>%
  ungroup() %>%
  group_by(species) %>%
  mutate(
    species = case_when(
      str_detect(species, "_pseudo_") ~ paste0(str_extract(species, "^[^_]+"), "_pseudo_species_", sprintf("%04d", cur_group_id())),
      TRUE ~ species
    )) %>%
  ungroup() %>%
  # Generate unique OTU_IDs and sum the abundances
  mutate(
    ASV_ID = OTU_ID,
    ASV_abundance = str_extract(OTU_ID, "(?<=size=)\\d+")
  ) %>%
  group_by(species) %>%
  mutate(
    OTU_abundance = sum(as.numeric(ASV_abundance))
  ) %>%
  arrange(desc(OTU_abundance)) %>%
  mutate(
    OTU_ID = paste0("OTU_", cur_group_id())
  ) %>%
  ungroup() %>%
  # Select reference_IDs for OTUs based on abundance
  group_by(OTU_ID, reference_ID) %>%
  mutate(
    reference_ID_abundance = sum(as.numeric(ASV_abundance))
  ) %>%
  ungroup() %>%
  arrange(desc(reference_ID_abundance)) %>%
  group_by(OTU_ID) %>%
  mutate(reference_ID = first(reference_ID)) %>%
  ungroup() %>%
  print(n = 500)

# Update the reference_IDs for closed-reference annotations ####################

# Pull the reference_IDs for closed-reference annotations
closed_reference_annotations <- taxonomy_1 %>%
  filter(
    str_detect(reference_ID, "ASV_")
  ) %>%
  unique(.) %>%
  pull(reference_ID)

# Pull representative OTUs
representative_references <- taxonomy_1 %>%
  filter(
    ASV_ID %in% closed_reference_annotations
  ) %>%
  select(rep_OTU_ID = OTU_ID, ASV_ID) %>%
  unique(.)

# Update the reference_IDs for closed-reference annotations
taxonomy_2 <- taxonomy_1 %>%
  # Join on ASV_ID
  left_join(representative_references, by = c("reference_ID" = "ASV_ID")) %>%
  mutate(
    reference_ID = case_when(
      # Replace ASV_ references with representative_OTU_ID if it exists
      str_detect(reference_ID, "ASV_") & !is.na(rep_OTU_ID) ~ rep_OTU_ID, 
      TRUE ~ reference_ID
    )
  ) %>%
  # Select the final columns
  select(
    OTU_ID, reference_ID, phylum, class, order, family, genus, species, score, cutoff, ASV_ID, ASV_abundance, OTU_abundance
  ) %>%
  print(n = 500)

# Generate an OTU table ########################################################

OTU_table_1 <- fread("data/prepared_reads/ASVs.txt") %>%
  select(ASV_ID = OTU_ID, everything()) %>%
  pivot_longer(
    cols = -ASV_ID,
    names_to = "sample_ID",
    values_to = "ASV_abundance"
  ) %>%
  left_join(
    taxonomy_2 %>%
      select(ASV_ID, OTU_ID) %>%
      mutate(ASV_ID = str_extract(ASV_ID, "^[^;]+")),
    by = "ASV_ID"
  ) %>%
  group_by(OTU_ID, sample_ID) %>%
  summarise(
    OTU_abundance = sum(ASV_abundance)
  ) %>%
  ungroup(.) %>%
  pivot_wider(
    names_from = sample_ID,
    values_from = OTU_abundance,
    values_fill = list(OTU_abundance = 0)
  ) %>%
  arrange(desc(rowSums(select(., -OTU_ID)))) %>%
  print(.)

# Sample-wise filtering: Filter low abundance OTUs sample-wise to control for
# wet lab contamination
source("code/abundance_filters.R")
OTU_table_2 <- filter_samples(OTU_table_1, 0.05)

# No change in the number of OTUs
nrow(OTU_table_1)
nrow(OTU_table_2)

# Library-wise filtering: Filter high abundance OTUs that occur in low abundance
# in specific samples to control for index switching
OTU_table_3 <- filter_library(OTU_table_2, 0.1)

# We lost 5836 or 0.1% of the total abundance or 0.24 % of the total reads
lost_reads <- sum(OTU_table_2 %>% select(-OTU_ID)) - sum(OTU_table_3 %>% select(-OTU_ID))
lost_reads / sum(OTU_table_2 %>% select(-OTU_ID)) * 100

# Compute OTU abundance
OTU_abundances <- OTU_table_3 %>%
  pivot_longer(
    cols = -OTU_ID,
    names_to = "sample_ID",
    values_to = "OTU_abundance"
  ) %>%
  group_by(OTU_ID) %>%
  summarise(
    abundance = sum(OTU_abundance)
  ) %>%
  print(.)

# Update the taxa table ########################################################

# Filter the taxonomy table to the filtered OTUs
taxonomy_3 <- taxonomy_2 %>%
  filter(OTU_ID %in% OTU_table_3$OTU_ID) %>%
  print(n = 500)

# Update the fasta file ########################################################

# Select the representative sequence for each OTU based on ASV abundance
fasta_headers <- taxonomy_2 %>%
  group_by(OTU_ID, phylum, class, order, family, genus, species) %>%
  arrange(desc(ASV_abundance)) %>%
  summarise(
    OTU_ID = first(OTU_ID),
    ASV_ID = first(ASV_ID)
  ) %>%
  ungroup(.) %>%
  mutate(
    taxonomy = str_c(phylum, class, order, family, genus, species, sep = ";"),
    header = paste0(OTU_ID, " ", taxonomy)
  ) %>%
  left_join(
    .,
    taxonomy_2 %>%
      select(OTU_ID, OTU_abundance),
    by = "OTU_ID"
  ) %>%
  unique(.) %>%
  arrange(desc(OTU_abundance)) %>%
  select(ASV_ID, header)
  
# Generate the final fasta file:
# Step 1: Read the FASTA file
fasta <- readDNAStringSet("data/prepared_reads/ASVs.fasta")

# Step 2: Filter the FASTA sequences by ASV_IDs in fasta_headers
filtered_fasta <- fasta[names(fasta) %in% fasta_headers$ASV_ID]

# Step 3: Rename the FASTA headers using the headers from fasta_headers
names(filtered_fasta) <- fasta_headers %>%
  filter(ASV_ID %in% names(fasta)) %>%
  pull(header)

# Save the final outputs #######################################################
fwrite(OTU_table_3, "data/prepared_reads/OTUs.txt", sep = "\t")
taxonomy_3 %>%
  filter(ASV_ID %in% fasta_headers[["ASV_ID"]]) %>%
  left_join(
    OTU_abundances %>%
      select(OTU_ID, abundance),
    by = "OTU_ID"
  ) %>%
  select(
    OTU_ID, reference_ID, phylum, class, order, family, genus, species, score, cutoff, abundance
  ) %>%
  unique(.) %>%
  arrange(desc(abundance)) %>%
  # Add fungal trait information
  left_join(
    .,
    fread("code/fungal_traits_genera.txt"),
    by = "genus"
  ) %>%
  fwrite(., "data/taxonomy/taxa_OTUs.txt", sep = "\t")
writeXStringSet(filtered_fasta, "data/prepared_reads/OTUs.fasta")
