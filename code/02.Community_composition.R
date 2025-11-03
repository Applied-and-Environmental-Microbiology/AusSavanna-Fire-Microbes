
# Required packages
require(vegan)
require(pairwiseAdonis)
require(tidyverse)

# Define the factor levels
litter_levels <- c("Leaf", "Grass")
regime_levels <- c("E2", "L2", "U")
community_levels <- c("Soil bacteria", "Litter bacteria" ,"Soil fungi", "Litter fungi")
regime_litter_type_levels <- c("E2_Leaf", "E2_Grass", "L2_Leaf", "L2_Grass", "U_Leaf", "U_Grass")

# Read in the sample metadata
metadata_fungi_soil <- data.table::fread("data/sample_metadata_fungi.txt") %>%
  filter(community == "Soil") %>%
  select(-c("community", "pco1_57_leaf_fungi", "pco2_38_leaf_fungi")) %>%
  # Order the regime and litter type factors
  mutate(
    regime = factor(regime, levels = regime_levels),
    litter_type = factor(litter_type, levels = litter_levels),
    regime_litter_type = factor(paste(regime, litter_type, sep = "_"), levels = regime_litter_type_levels)
  )
metadata_fungi_litter <- data.table::fread("data/sample_metadata_fungi.txt") %>%
  filter(community == "Leaf") %>%
  select(-c("community", "pco1_54_soil_fungi", "pco2_39_soil_fungi")) %>%
  # Order the regime and litter type factors
  mutate(
    regime = factor(regime, levels = regime_levels),
    litter_type = factor(litter_type, levels = litter_levels),
    regime_litter_type = factor(paste(regime, litter_type, sep = "_"), levels = regime_litter_type_levels)
  )
metadata_bacteria_soil <- data.table::fread("data/sample_metadata_bacteria.txt") %>%
  filter(community == "Soil") %>%
  select(-c("community", "pco1_49_leaf_bacteria", "pco2_42_leaf_bacteria")) %>%
  # Order the regime and litter type factors
  mutate(
    regime = factor(regime, levels = regime_levels),
    litter_type = factor(litter_type, levels = litter_levels),
    regime_litter_type = factor(paste(regime, litter_type, sep = "_"), levels = regime_litter_type_levels)
  )
metadata_bacteria_litter <- data.table::fread("data/sample_metadata_bacteria.txt") %>%
  filter(community == "Leaf") %>%
  select(-c("community", "pco1_54_soil_bacteria", "pco2_40_soil_bacteria")) %>%
  # Order the regime and litter type factors
  mutate(
    regime = factor(regime, levels = regime_levels),
    litter_type = factor(litter_type, levels = litter_levels),
    regime_litter_type = factor(paste(regime, litter_type, sep = "_"), levels = regime_litter_type_levels)
  )

# Read in the distance matices
dist_fungi_soil <- as.dist(
  data.table::fread("data/distance_fungi_soil.txt") %>%
    column_to_rownames("V1")
  )
dist_fungi_litter <- as.dist(
  data.table::fread("data/distance_fungi_leaf.txt") %>%
    column_to_rownames("V1")
  )
dist_bacteria_soil <- as.dist(
  data.table::fread("data/distance_bacteria_soil.txt") %>%
    column_to_rownames("V1")
  )
dist_bacteria_litter <- as.dist(
  data.table::fread("data/distance_bacteria_leaf.txt") %>%
    column_to_rownames("V1")
  )

# (1) Run PERMANOVA ------------------------------------------------------------

# Set seed for reproducibility
set.seed(1986)

#### (1a) PERMANOVAs ####

# Loop through each combination of community and environment

# Initialise a list
permanova_results <- list()

# Loop through each combination
for (community in c("fungi", "bacteria")) {
  for (environment in c("soil", "litter")) {
    
    # Select the appropriate distance matrix and metadata
    dist_matrix <- get(paste0("dist_", community, "_", environment))
    metadata <- get(paste0("metadata_", community, "_", environment))
    
    # Run PERMANOVA
    permanova <- adonis2(
      dist_matrix ~ regime * litter_type,
      data = metadata,
      by = "terms",
      strata = metadata$block
    )
    
    # Store the results in the list
    permanova_results[[paste0(community, "_", environment)]] <- permanova
    
    # Pairwise PERMANOVA for different regimes
    pairwise_regime <- pairwise.adonis(dist_matrix, metadata$regime, p.adjust.m = "none")
    permanova_results[[paste0(community, "_", environment, "_regime")]] <- pairwise_regime
    
    # Pairwise PERMANOVA for different litter types
    pairwise_litter <- pairwise.adonis(dist_matrix, metadata$litter_type, p.adjust.m = "none")
    permanova_results[[paste0(community, "_", environment, "_litter")]] <- pairwise_litter
    
    # Pairwise PERMANOVA for different regime and litter type combinations
    pairwise_regime_litter <- pairwise.adonis(dist_matrix, metadata$regime_litter_type, p.adjust.m = "none")
    permanova_results[[paste0(community, "_", environment, "_regime_litter")]] <- pairwise_regime_litter
    
  }
}

#### (1b) Format the results ####

# Create a table of main PERMANOVA results
main_permanova_table <- bind_rows(
  lapply(c("fungi_soil", "fungi_litter", "bacteria_soil", "bacteria_litter"), function(name) {
    result <- permanova_results[[name]]
    result <- as.data.frame(result)
    result$Factor <- rownames(result)
    result$Community <- name
    rownames(result) <- NULL
    return(result)
  })
) %>%
  select(
    Community, Factor, `pseudo-F` = `F`, `R2 (%)` = R2, `P-value` = `Pr(>F)`
  ) %>%
  mutate(
    `pseudo-F` = formatC(`pseudo-F`, format = "f", digits = 2),
    `R2 (%)` = formatC(`R2 (%)` * 100, format = "f", digits = 0),
    `P-value` = formatC(`P-value`, format = "f", digits = 3),
    Community = case_when(
      Community == "fungi_soil" ~ "Soil fungi",
      Community == "fungi_litter" ~ "Litter fungi",
      Community == "bacteria_soil" ~ "Soil bacteria",
      Community == "bacteria_litter" ~ "Litter bacteria",
      TRUE ~ Community
    ),
    Factor = case_when(
      Factor == "regime" ~ "Fire regime",
      Factor == "litter_type" ~ "Litter type",
      Factor == "regime:litter_type" ~ "Interaction",
      TRUE ~ Factor
    )) %>%
  filter(!Factor %in% c("Total", "Residual")) %>%
  # Level the community factor
  mutate(Community = factor(Community, levels = community_levels)) %>%
  arrange(Community)

# Create table for pairwise comparisons of fire regimes
pairwise_regime_table <- bind_rows(
  lapply(c("fungi_soil_regime", "fungi_litter_regime", "bacteria_soil_regime", "bacteria_litter_regime"), function(name) {
    result <- permanova_results[[name]]
    result <- as.data.frame(result)
    result$Comparison <- rownames(result)
    result$Community <- name
    rownames(result) <- NULL
    return(result)
  })
) %>%
  select(
    Community, Comparison = pairs, `pseudo-F` = F.Model, `R2 (%)` = R2, `P-value` = `p.adjusted`
  ) %>%
  mutate(
    `pseudo-F` = formatC(`pseudo-F`, format = "f", digits = 2),
    `R2 (%)` = formatC(`R2 (%)` * 100, format = "f", digits = 0),
    `P-value` = formatC(`P-value`, format = "f", digits = 3),
    Community = case_when(
      Community == "fungi_soil_regime" ~ "Soil fungi",
      Community == "fungi_litter_regime" ~ "Litter fungi",
      Community == "bacteria_soil_regime" ~ "Soil bacteria",
      Community == "bacteria_litter_regime" ~ "Litter bacteria",
      TRUE ~ Community
    ),
    Comparison = case_when(
      Comparison == "E2 vs L2" ~ "Early vs. late",
      Comparison == "L2 vs E2" ~ "Early vs. late",
      Comparison == "E2 vs U" ~ "Early vs. unburnt",
      Comparison == "U vs E2" ~ "Early vs. unburnt",
      Comparison == "L2 vs U" ~ "Late vs. unburnt",
      Comparison == "U vs L2" ~ "Late vs. unburnt",
      TRUE ~ Comparison
    )
  ) %>%
  # Level the community and comparison factors
  mutate(
    Community = factor(Community, levels = community_levels),
    Comparison = factor(Comparison, levels = c("Early vs. late", "Early vs. unburnt", "Late vs. unburnt"))
    ) %>%
  arrange(Comparison) %>%
  arrange(Community)

# Create table for pairwise comparisons of fire regime-litter type combinations
pairwise_regime_litter_table <- bind_rows(
  lapply(c("fungi_soil_regime_litter", "fungi_litter_regime_litter", 
           "bacteria_soil_regime_litter", "bacteria_litter_regime_litter"), function(name) {
             result <- permanova_results[[name]]
             result <- as.data.frame(result)
             result$Comparison <- rownames(result)
             result$Community <- name
             rownames(result) <- NULL
             return(result)
           })
) %>%
  select(
    Community, Comparison = pairs, `pseudo-F` = F.Model, `R2 (%)` = R2, `P-value` = `p.adjusted`
  ) %>%
  mutate(
    `pseudo-F` = formatC(`pseudo-F`, format = "f", digits = 2),
    `R2 (%)`   = formatC(`R2 (%)` * 100, format = "f", digits = 0),
    `P-value`  = formatC(`P-value`, format = "f", digits = 3),
    Community  = case_when(
      Community == "fungi_soil_regime_litter"    ~ "Soil fungi",
      Community == "fungi_litter_regime_litter"  ~ "Litter fungi",
      Community == "bacteria_soil_regime_litter" ~ "Soil bacteria",
      Community == "bacteria_litter_regime_litter" ~ "Litter bacteria",
      TRUE ~ Community
    ),
    Comparison = case_when(
      # --- Early vs Late (all permutations) ---
      Comparison == "E2_Leaf vs L2_Leaf"   ~ "Early (leaf) vs. late (leaf)",
      Comparison == "L2_Leaf vs E2_Leaf"   ~ "Early (leaf) vs. late (leaf)",
      Comparison == "E2_Leaf vs L2_Grass"  ~ "Early (leaf) vs. late (grass)",
      Comparison == "L2_Grass vs E2_Leaf"  ~ "Early (leaf) vs. late (grass)",
      Comparison == "E2_Grass vs L2_Leaf"  ~ "Early (grass) vs. late (leaf)",
      Comparison == "L2_Leaf vs E2_Grass"  ~ "Early (grass) vs. late (leaf)",
      Comparison == "E2_Grass vs L2_Grass" ~ "Early (grass) vs. late (grass)",
      Comparison == "L2_Grass vs E2_Grass" ~ "Early (grass) vs. late (grass)",
      
      # --- Early vs Unburnt ---
      Comparison == "E2_Leaf vs U_Leaf"   ~ "Early (leaf) vs. unburnt (leaf)",
      Comparison == "U_Leaf vs E2_Leaf"   ~ "Early (leaf) vs. unburnt (leaf)",
      Comparison == "E2_Leaf vs U_Grass"  ~ "Early (leaf) vs. unburnt (grass)",
      Comparison == "U_Grass vs E2_Leaf"  ~ "Early (leaf) vs. unburnt (grass)",
      Comparison == "E2_Grass vs U_Leaf"  ~ "Early (grass) vs. unburnt (leaf)",
      Comparison == "U_Leaf vs E2_Grass"  ~ "Early (grass) vs. unburnt (leaf)",
      Comparison == "E2_Grass vs U_Grass" ~ "Early (grass) vs. unburnt (grass)",
      Comparison == "U_Grass vs E2_Grass" ~ "Early (grass) vs. unburnt (grass)",
      
      # --- Late vs Unburnt ---
      Comparison == "L2_Leaf vs U_Leaf"   ~ "Late (leaf) vs. unburnt (leaf)",
      Comparison == "U_Leaf vs L2_Leaf"   ~ "Late (leaf) vs. unburnt (leaf)",
      Comparison == "L2_Leaf vs U_Grass"  ~ "Late (leaf) vs. unburnt (grass)",
      Comparison == "U_Grass vs L2_Leaf"  ~ "Late (leaf) vs. unburnt (grass)",
      Comparison == "L2_Grass vs U_Leaf"  ~ "Late (grass) vs. unburnt (leaf)",
      Comparison == "U_Leaf vs L2_Grass"  ~ "Late (grass) vs. unburnt (leaf)",
      Comparison == "L2_Grass vs U_Grass" ~ "Late (grass) vs. unburnt (grass)",
      Comparison == "U_Grass vs L2_Grass" ~ "Late (grass) vs. unburnt (grass)",
      
      # --- Leaf vs Grass within same fire regime (always leaf first) ---
      Comparison == "E2_Leaf vs E2_Grass" ~ "Early (leaf) vs. early (grass)",
      Comparison == "E2_Grass vs E2_Leaf" ~ "Early (leaf) vs. early (grass)",
      Comparison == "L2_Leaf vs L2_Grass" ~ "Late (leaf) vs. late (grass)",
      Comparison == "L2_Grass vs L2_Leaf" ~ "Late (leaf) vs. late (grass)",
      Comparison == "U_Leaf vs U_Grass"   ~ "Unburnt (leaf) vs. unburnt (grass)",
      Comparison == "U_Grass vs U_Leaf"   ~ "Unburnt (leaf) vs. unburnt (grass)",
      
      TRUE ~ Comparison
    )
  ) %>%
  # Order factors
  mutate(
    Community = factor(Community, levels = community_levels),
    Comparison = factor(Comparison, levels = c(
      "Early (leaf) vs. late (leaf)",
      "Early (leaf) vs. late (grass)",
      "Early (grass) vs. late (leaf)",
      "Early (grass) vs. late (grass)",
      "Early (leaf) vs. unburnt (leaf)",
      "Early (leaf) vs. unburnt (grass)",
      "Early (grass) vs. unburnt (leaf)",
      "Early (grass) vs. unburnt (grass)",
      "Late (leaf) vs. unburnt (leaf)",
      "Late (leaf) vs. unburnt (grass)",
      "Late (grass) vs. unburnt (leaf)",
      "Late (grass) vs. unburnt (grass)",
      "Early (leaf) vs. early (grass)",
      "Late (leaf) vs. late (grass)",
      "Unburnt (leaf) vs. unburnt (grass)"
    ))
  ) %>%
  arrange(Comparison) %>%
  arrange(Community)

# View the results
main_permanova_table
pairwise_regime_table
pairwise_regime_litter_table

#### (1c) Save the results ####

data.table::fwrite(
  main_permanova_table,
  "output/permanova_main_results.txt",
  sep = "\t"
)
data.table::fwrite(
  pairwise_regime_table,
  "output/permanova_pairwise_results.txt",
  sep = "\t"
)
data.table::fwrite(
  pairwise_regime_litter_table,
  "output/permanova_pairwise_regime_litter_results.txt",
  sep = "\t"
)
