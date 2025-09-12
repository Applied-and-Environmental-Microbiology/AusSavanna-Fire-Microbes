
# Load packages
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(tidyverse))

# Step 1: Species and genus annotations
species <- fread(paste0(taxa_dir, "species_classification.txt")) %>%
  filter(genus != "unidentified")
species_ID <- species$ID

# Step 2: Genus and family annotations
genus <- fread(paste0(taxa_dir, "genus_classification.txt")) %>%
  filter(
    family != "unidentified",
    !ID %in% species_ID
  )
genus_ID <- genus$ID
species_2 <- fread(paste0(taxa_dir, "species_classification.txt")) %>%
  filter(
    family != "unidentified",
    !ID %in% c(species_ID, genus_ID)
  )
species_ID_2 <- species_2$ID

# Sterp 3: Family and order annotations
family <- fread(paste0(taxa_dir, "family_classification.txt")) %>%
  filter(
    order != "unidentified",
    !ID %in% c(species_ID, genus_ID, species_ID_2)
  )
family_ID <- family$ID
species_3 <- fread(paste0(taxa_dir, "species_classification.txt")) %>%
  filter(
    order != "unidentified",
    !ID %in% c(species_ID, genus_ID, species_ID_2, family_ID)
  )
species_ID_3 <- species_3$ID
genus_2 <- fread(paste0(taxa_dir, "genus_classification.txt")) %>%
  filter(
    order != "unidentified",
    !ID %in% c(species_ID, genus_ID, species_ID_2, family_ID, species_ID_3)
  )
genus_ID_2 <- genus_2$ID

# Step 4: Order and class annotations
order <- fread(paste0(taxa_dir, "order_classification.txt")) %>%
  filter(
    class != "unidentified",
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2
      )
  )
order_ID <- order$ID
species_4 <- fread(paste0(taxa_dir, "species_classification.txt")) %>%
  filter(
    class != "unidentified",
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2,
      order_ID
      )
  )
species_ID_4 <- species_4$ID
genus_3 <- fread(paste0(taxa_dir, "genus_classification.txt")) %>%
  filter(
    class != "unidentified",
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2,
      order_ID, species_ID_4
      )
  )
genus_ID_3 <- genus_3$ID
family_2 <- fread(paste0(taxa_dir, "family_classification.txt")) %>%
  filter(
    class != "unidentified",
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2,
      order_ID, species_ID_4, genus_ID_3
      )
  )
family_ID_2 <- family_2$ID

# Sterp 5: Class and phylum annotations
class <- fread(paste0(taxa_dir, "class_classification.txt")) %>%
  filter(
    phylum != "unidentified",
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2,
      order_ID, species_ID_4, genus_ID_3, family_ID_2
      )
  )
class_ID <- class$ID
species_5 <- fread(paste0(taxa_dir, "species_classification.txt")) %>%
  filter(
    phylum != "unidentified",
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2,
      order_ID, species_ID_4, genus_ID_3, family_ID_2, class_ID
      )
  )
species_ID_5 <- species_5$ID
genus_4 <- fread(paste0(taxa_dir, "genus_classification.txt")) %>%
  filter(
    phylum != "unidentified",
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2,
      order_ID, species_ID_4, genus_ID_3, family_ID_2, class_ID, species_ID_5
      )
  )
genus_ID_4 <- genus_4$ID
family_3 <- fread(paste0(taxa_dir, "family_classification.txt")) %>%
  filter(
    phylum != "unidentified",
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2,
      order_ID, species_ID_4, genus_ID_3, family_ID_2, class_ID, species_ID_5,
      genus_ID_4
      )
  )
family_ID_3 <- family_3$ID
order_2 <- fread(paste0(taxa_dir, "order_classification.txt")) %>%
  filter(
    phylum != "unidentified",
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2,
      order_ID, species_ID_4, genus_ID_3, family_ID_2, class_ID, species_ID_5,
      genus_ID_4, family_ID_3
      )
  )
order_ID_2 <- order_2$ID

# Step 6: Phylum and kingdom annotations
all <- fread(paste0(taxa_dir, "phylum_classification.txt")) %>%
  filter(
    !ID %in% c(
      species_ID, genus_ID, species_ID_2, family_ID, species_ID_3, genus_ID_2,
      order_ID, species_ID_4, genus_ID_3, family_ID_2, class_ID, species_ID_5,
      genus_ID_4, family_ID_3, order_ID_2
    )
  )

# Build the taxonomy table
taxonomy <- rbind(
  species, genus, species_2, family, species_3, genus_2, order, species_4,
  genus_3, family_2, class, species_5, genus_4, family_3, order_2, all
) %>%
  mutate(
    OTU_ID = ID,
    confidence = as.numeric(confidence),
    cutoff = as.numeric(cutoff),
    score = as.numeric(score)
  ) %>%
  select(OTU_ID, everything(), -ID) %>%
  arrange(
    desc(confidence),
    desc(score)
  )

# Check if any duplicate OTU_IDs
duplicates <- taxonomy %>%
  group_by(OTU_ID) %>%
  filter(n() > 1) %>%
  ungroup()

# Print duplicates
print(duplicates, n = Inf)

# Exit with an error if duplicates are found
if (nrow(duplicates) > 0) {
  stop("Error: Duplicate OTU_IDs found.")
}

# Save the taxonomy table
taxonomy %>%
  # Replace spaces with underscores in all columns
  mutate_all(~ gsub(" ", "_", .)) %>%
  fwrite(., paste0(taxa_dir, "/taxa_ASVs.txt"), sep = "\t")