# Removes OTUs sample-wise with an abundance <= threshold

filter_samples <- function(otu_tab, threshold) {
  # Filter out zero values and calculate relative abundance
  otu_tab %>%
    pivot_longer(-OTU_ID, names_to = 'sample') %>%
    filter(value > 0) %>%
    group_by(sample) %>%
    mutate(rel_abund = value / sum(value) * 100) %>%
    ungroup() %>%
    
    # Set values below the threshold to zero
    mutate(value = replace(value, rel_abund <= threshold, 0)) %>%
    
    # Remove rows where all values are zero 
    select(-rel_abund) %>%
    filter(value > 0) %>%
    
    # Reshape the table
    pivot_wider(names_from = "sample", values_from = "value", values_fill = 0) #%>%
    #column_to_rownames(var = "OTU_ID")
  
}

#### Filter rare occurrences of abundant OTUs ###################################

# Proportional version
filter_library <- function(library_specific_OTU_table, threshold) {
  
  # Remove rare occurrences of abundant OTUs using a 0.5% threshold.
  filtered_otu_table <- library_specific_OTU_table %>%
    #rownames_to_column(var = "OTU_ID") %>%
    pivot_longer(-OTU_ID, names_to = "sample",
                 values_to = "OTU_sample_abundance") %>%
    group_by(OTU_ID) %>%
    mutate(
      OTU_library_abundance = sum(OTU_sample_abundance),
      rel_abundance = (OTU_sample_abundance / OTU_library_abundance) * 100
    ) %>%
    mutate(OTU_sample_abundance = replace(
      OTU_sample_abundance, rel_abundance <= threshold, 0)) %>%
    ungroup() %>%
    select(-c(OTU_library_abundance, rel_abundance)) %>%
    pivot_wider(
      names_from = "sample",
      values_from = "OTU_sample_abundance",
      values_fill = 0
      )
  
  return(filtered_otu_table)
  
}