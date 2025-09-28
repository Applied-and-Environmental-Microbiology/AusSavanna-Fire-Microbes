
require(data.table)
require(tidyverse)
source("code/abundance_filters.R")

ASVs <- fread("data/ASVs.txt") %>%
  rename(OTU_ID = ASV_ID)
taxa <- fread("data/taxonomy.txt")

# Filter ASVs
ASVs_filtered_1 <- filter_samples(ASVs, 0.05)
ASVs_filtered_2 <- filter_library(ASVs_filtered_1, 0.5)

nrow(ASVs)
nrow(ASVs_filtered_1)

sum(ASVs %>% select(-OTU_ID))
sum(ASVs_filtered_1 %>% select(-OTU_ID))
sum(ASVs_filtered_2 %>% select(-OTU_ID))

# Filter taxa to ASVs
taxa_filtered <- taxa %>%
  filter(ASV_ID %in% ASVs_filtered_2$OTU_ID)

# Save filtered ASVs and taxa
ASVs_filtered_2 %>%
  rename(ASV_ID = OTU_ID) %>%
  fwrite(., "data/ASVs_abundance_filtered.txt")
fwrite(taxa_filtered, "data/taxonomy_abundance_filtered.txt")
