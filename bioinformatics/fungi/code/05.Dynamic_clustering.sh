#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=8G
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=slurm/%x.%j.out

# Activate the conda environment
echo "Starting at:" $(date)
source ~/.bashrc
# Activate the conda environment
conda activate dynamic_clustering

Rscript --vanilla - <<EOF
# Required packages
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(tidyverse))

# Read in the data
taxa_file <- fread("../data/taxonomy/taxa_ASVs.txt") %>%
  dplyr::rename(reference_ID = ReferenceID) %>%
  print(.)
cutoff_file <- fread("./dnabarcoder/unite2024ITS2.unique.cutoffs.best.txt") %>%
  select(rank = Rank, taxa = Dataset, cutoff = "cut-off") %>%
  print(.)
asv_sequences <- readDNAStringSet("../data/prepared_reads/ASVs.fasta") %>%
  print(.)

# Constants
threads <- 24 # Number of threads
minlen <- 50  # Minimum BLAST alignment length
files <- list.files("tmp", full.names = TRUE) # Files to clear after each rank

# Create temporary directories for intermediate files
dir.create("tmp")
dir.create("tmp_clusters")

# Match unique taxa to cut-offs
unique_taxa_cutoffs <- taxa_file %>%
  select(-c("OTU_ID", "reference_ID", "score", "cutoff", "confidence")) %>%
  pivot_longer(
    -"rank", 
    names_to = "level",
    values_to = "taxa"
  ) %>%
  filter(
    !taxa %in% c("unidentified", "") & !is.na(taxa)
  ) %>%
  select(-level) %>%
  unique(.) %>%
  left_join(
    .,
    cutoff_file,
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  print(.)

# Define ranks and superranks for looping
ranks <- c("species", "genus", "family", "order", "class", "phylum")
superranks_species <- c("genus", "family", "order", "class", "phylum")
superranks_genus <- c("family", "order", "class", "phylum")
superranks_family <- c("order", "class", "phylum")
superranks_order <- c("class", "phylum")
superranks_class <- c("phylum")

# !!! Do better here !!!
# Make these pseudo cutoffs taxanomically infromed next
pseudo_cutoffs <- cutoff_file %>%
  filter(
    rank == "species" & taxa == "All" |
      rank == "genus" & taxa == "All" |
      rank == "family" & taxa == "All" |
      rank == "order" & taxa == "All" |
      rank == "clasee" & taxa == "All" |
      rank == "phylum" & taxa == "All"
  )

###############################################################################*
# (1) Match cut-offs to classified ASVs #######################################
###############################################################################*

### (1.1) Species ##############################################################

# Match cutoffs to unique taxa for each rank
species_cutoffs_all <- taxa_file %>%
  select(
    rank, all_of(superranks_species)
  ) %>%
  filter(rank == "species") %>%
  pivot_longer(
    -rank,
    values_to = "taxa"
  ) %>%
  left_join(
    .,
    unique_taxa_cutoffs %>%
      filter(rank == "species"),
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  select(-rank, rank = name) %>%
  unique(.)

# Create a list to store the filtered data frames
species_cutoffs <- list()

for (superrank in superranks_species) {
  species_cutoffs[[superrank]] <- species_cutoffs_all %>% 
    filter(rank == superrank) %>%
    select(-rank)
}

# Assign species cutoffs to the taxa file
taxa_file_1 <- taxa_file %>%
  mutate(
    species_cutoff = case_when(
      genus %in% species_cutoffs[["genus"]][["taxa"]] ~ species_cutoffs[["genus"]][["cutoff"]][match(genus, species_cutoffs[["genus"]][["taxa"]])],
      family %in% species_cutoffs[["family"]][["taxa"]] ~ species_cutoffs[["family"]][["cutoff"]][match(family, species_cutoffs[["family"]][["taxa"]])],
      order %in% species_cutoffs[["order"]][["taxa"]] ~ species_cutoffs[["order"]][["cutoff"]][match(order, species_cutoffs[["order"]][["taxa"]])],
      class %in% species_cutoffs[["class"]][["taxa"]] ~ species_cutoffs[["class"]][["cutoff"]][match(class, species_cutoffs[["class"]][["taxa"]])],
      phylum %in% species_cutoffs[["phylum"]][["taxa"]] ~ species_cutoffs[["phylum"]][["cutoff"]][match(phylum, species_cutoffs[["phylum"]][["taxa"]])],
      TRUE ~ cutoff_file %>% filter(rank == "species" & taxa == "All") %>% pull(cutoff)
    )
  )

### (1.2) Genus ################################################################

# Match cutoffs to unique taxa for each rank
genus_cutoffs_all <- taxa_file %>%
  select(
    rank, all_of(superranks_genus)
  ) %>%
  filter(rank == "genus") %>%
  pivot_longer(
    -rank,
    values_to = "taxa"
  ) %>%
  left_join(
    .,
    unique_taxa_cutoffs %>%
      filter(rank == "genus"),
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  select(-rank, rank = name) %>%
  unique(.)

# Create a list to store the filtered data frames
genus_cutoffs <- list()

for (superrank in superranks_genus) {
  genus_cutoffs[[superrank]] <- genus_cutoffs_all %>% 
    filter(rank == superrank) %>%
    select(-rank)
}

# Assign genus cutoffs to the taxa file
taxa_file_2 <- taxa_file_1 %>%
  mutate(
    genus_cutoff = case_when(
      family %in% genus_cutoffs[["family"]][["taxa"]] ~ genus_cutoffs[["family"]][["cutoff"]][match(family, genus_cutoffs[["family"]][["taxa"]])],
      order %in% genus_cutoffs[["order"]][["taxa"]] ~ genus_cutoffs[["order"]][["cutoff"]][match(order, genus_cutoffs[["order"]][["taxa"]])],
      class %in% genus_cutoffs[["class"]][["taxa"]] ~ genus_cutoffs[["class"]][["cutoff"]][match(class, genus_cutoffs[["class"]][["taxa"]])],
      phylum %in% genus_cutoffs[["phylum"]][["taxa"]] ~ genus_cutoffs[["phylum"]][["cutoff"]][match(phylum, genus_cutoffs[["phylum"]][["taxa"]])],
      TRUE ~ cutoff_file %>% filter(rank == "genus" & taxa == "All") %>% pull(cutoff)
    )
  )

### (1.3) Family ###############################################################

# Match cutoffs to unique taxa for each rank
family_cutoffs_all <- taxa_file %>%
  select(
    rank, all_of(superranks_family)
  ) %>%
  filter(rank == "family") %>%
  pivot_longer(
    -rank,
    values_to = "taxa"
  ) %>%
  left_join(
    .,
    unique_taxa_cutoffs %>%
      filter(rank == "family"),
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  select(-rank, rank = name) %>%
  unique(.)

# Create a list to store the filtered data frames
family_cutoffs <- list()

for (superrank in superranks_family) {
  family_cutoffs[[superrank]] <- family_cutoffs_all %>% 
    filter(rank == superrank) %>%
    select(-rank)
}

# Assign family cutoffs to the taxa file
taxa_file_3 <- taxa_file_2 %>%
  mutate(
    family_cutoff = case_when(
      order %in% family_cutoffs[["order"]][["taxa"]] ~ family_cutoffs[["order"]][["cutoff"]][match(order, family_cutoffs[["order"]][["taxa"]])],
      class %in% family_cutoffs[["class"]][["taxa"]] ~ family_cutoffs[["class"]][["cutoff"]][match(class, family_cutoffs[["class"]][["taxa"]])],
      phylum %in% family_cutoffs[["phylum"]][["taxa"]] ~ family_cutoffs[["phylum"]][["cutoff"]][match(phylum, family_cutoffs[["phylum"]][["taxa"]])],
      TRUE ~ cutoff_file %>% filter(rank == "family" & taxa == "All") %>% pull(cutoff)
    )
  )

### (1.4) Order ################################################################

# Match cutoffs to unique taxa for each rank
order_cutoffs_all <- taxa_file %>%
  select(
    rank, all_of(superranks_order)
  ) %>%
  filter(rank == "order") %>%
  pivot_longer(
    -rank,
    values_to = "taxa"
  ) %>%
  left_join(
    .,
    unique_taxa_cutoffs %>%
      filter(rank == "order"),
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  select(-rank, rank = name) %>%
  unique(.)

# Create a list to store the filtered data frames
order_cutoffs <- list()

for (superrank in superranks_order) {
  order_cutoffs[[superrank]] <- order_cutoffs_all %>% 
    filter(rank == superrank) %>%
    select(-rank)
}

# Assign order cutoffs to the taxa file
taxa_file_4 <- taxa_file_3 %>%
  mutate(
    order_cutoff = case_when(
      class %in% order_cutoffs[["class"]][["taxa"]] ~ order_cutoffs[["class"]][["cutoff"]][match(class, order_cutoffs[["class"]][["taxa"]])],
      phylum %in% order_cutoffs[["phylum"]][["taxa"]] ~ order_cutoffs[["phylum"]][["cutoff"]][match(phylum, order_cutoffs[["phylum"]][["taxa"]])],
      TRUE ~ cutoff_file %>% filter(rank == "order" & taxa == "All") %>% pull(cutoff)
    )
  )

### (1.5) Class ################################################################

# Match cutoffs to unique taxa for each rank
class_cutoffs_all <- taxa_file %>%
  select(
    rank, all_of(superranks_class)
  ) %>%
  filter(rank == "class") %>%
  pivot_longer(
    -rank,
    values_to = "taxa"
  ) %>%
  left_join(
    .,
    unique_taxa_cutoffs %>%
      filter(rank == "class"),
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  select(-rank, rank = name) %>%
  unique(.)

# Create a list to store the filtered data frames
class_cutoffs <- list()

for (superrank in superranks_class) {
  class_cutoffs[[superrank]] <- class_cutoffs_all %>% 
    filter(rank == superrank) %>%
    select(-rank)
}

# Assign class cutoffs to the taxa file
taxa_file_5 <- taxa_file_4 %>%
  mutate(
    class_cutoff = case_when(
      phylum %in% class_cutoffs[["phylum"]][["taxa"]] ~ class_cutoffs[["phylum"]][["cutoff"]][match(phylum, class_cutoffs[["phylum"]][["taxa"]])],
      TRUE ~ cutoff_file %>% filter(rank == "class" & taxa == "All") %>% pull(cutoff)
    )
  )

### (1.6) Phylum ###############################################################

# Match cutoffs to unique taxa for each rank
phylum_cutoffs_all <- taxa_file %>%
  select(
    rank, kingdom
  ) %>%
  filter(rank == "phylum") %>%
  pivot_longer(
    -rank,
    values_to = "taxa"
  ) %>%
  left_join(
    .,
    unique_taxa_cutoffs %>%
      filter(rank == "phylum"),
    by = c("rank", "taxa")
  ) %>%
  filter(
    !is.na(cutoff)
  ) %>%
  select(-rank, rank = name) %>%
  unique(.)

# Create a list to store the filtered data frames
phylum_cutoffs <- list()

# No superranks for phylum, so just store the data frame
phylum_cutoffs[["phylum"]] <- phylum_cutoffs_all %>% 
  select(-rank)

# Assign phylum cutoffs to the taxa file
taxa_cutoffs <- taxa_file_5 %>%
  mutate(
    phylum_cutoff = case_when(
      phylum %in% phylum_cutoffs[["phylum"]][["taxa"]] ~ phylum_cutoffs[["phylum"]][["cutoff"]][match(phylum, phylum_cutoffs[["phylum"]][["taxa"]])],
      TRUE ~ cutoff_file %>% filter(rank == "phylum" & taxa == "All") %>% pull(cutoff)
    )
  )

### (1.7) Clean-up intermediate files and variables ############################

suppressWarnings({
  for (rank in ranks) {
    rm(list = ls(pattern = (paste0(rank, "_cutoff"))))
  }

  rm(list = c(
    ls(pattern = "taxa_file"),
    ls(pattern = "superrank"),
    "rank",
    "unique_taxa_cutoffs",
    "cutoff_file"
  ))
})

###############################################################################*
# (2) Phylum clusters #########################################################
###############################################################################*

# Define the ranks of interest and its enclosing supertaxon
this_superrank <- "kingdom"
this_rank <- "phylum"
this_subrank <- "class"

# Write the taxonomy file as the phylum file
taxa_cutoffs %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Kingdom cut-offs
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (2a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]
  
  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))
  
  # Check if this_cutoff is 0 or NA (indicating a pseudo taxa)
  if (is.na(this_cutoff) || this_cutoff == 0) {
    # Get cutoff for pseudo taxa from pseudo_cutoffs
    this_cutoff <- pseudo_cutoffs %>%
      filter(rank == this_rank) %>%
      pull(cutoff)
    
    message(paste0("Assigning cutoff ", this_cutoff, " to pseudo taxa ", this_supertaxon))
  }
  
  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
    
    identified_asvs <- taxa_cutoffs %>%
      #filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }
    
    unidentified_asvs <- taxa_cutoffs %>%
      #filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }
    
    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]
    
    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")
    
    ## Cluster unidentified ASVs ####
    
    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))
    
    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )
    
    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length bitscore qstart qend'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )
    
    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        coverage = V4,
        bitscore = V5,
        start = V6,
        end = V7
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score with the necessary adjustment
      mutate(
        score = sim / 100,
        score = case_when(
          coverage < minlen ~ (score * coverage) / minlen,
          TRUE ~ score
        )
      ) %>%
      filter(score > this_cutoff) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        kingdom = this_supertaxon,
        class = "unidentified",
        order = "unidentified",
        family = "unidentified",
        genus = "unidentified",
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        confidence = 0,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        kingdom, phylum, class, order, family, genus, species,
        rank, cutoff, confidence, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      )
    
    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )
    
    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)
    
    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
    
    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }
  
  ### (3b) De novo clustering ###########################################
  
  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)
  
  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))
    
    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")
    
    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))
    
    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")
    
    # Initialise pseudo_id
    pseudo_id <- 1
    
    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)
    
    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (3) Class clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "phylum"
this_rank <- "class"
this_subrank <- "order"

# Write the class file as the order file
fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Class cut-offs
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (3a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]
  
  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))
  
  # Check if this_cutoff is 0 or NA (indicating a pseudo taxa)
  if (is.na(this_cutoff) || this_cutoff == 0) {
    # Get cutoff for pseudo taxa from pseudo_cutoffs
    this_cutoff <- pseudo_cutoffs %>%
      filter(rank == this_rank) %>%
      pull(cutoff)
    
    message(paste0("Assigning cutoff ", this_cutoff, " to pseudo taxa ", this_supertaxon))
  }
  
  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
    
    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }
    
    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }
    
    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]
    
    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")
    
    ## Cluster unidentified ASVs ####
    
    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))
    
    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )
    
    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length bitscore qstart qend'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )
    
    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        coverage = V4,
        bitscore = V5,
        start = V6,
        end = V7
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score with the necessary adjustment
      mutate(
        score = sim / 100,
        score = case_when(
          coverage < minlen ~ (score * coverage) / minlen,
          TRUE ~ score
        )
      ) %>%
      filter(score > this_cutoff) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        phylum = this_supertaxon,
        order = "unidentified",
        family = "unidentified",
        genus = "unidentified",
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        confidence = 0,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        phylum, class, order, family, genus, species,
        rank, cutoff, confidence, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      )
    
    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )
    
    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)
    
    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
    
    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }
  
  ### (3b) De novo clustering ###########################################
  
  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)
  
  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))
    
    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")
    
    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))
    
    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")
    
    # Initialise pseudo_id
    pseudo_id <- 1
    
    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)
    
    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (4) Order clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "class"
this_rank <- "order"
this_subrank <- "family"

# Write the class file as the order file
fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Class cut-offs
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (4a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]
  
  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))
  
  # Check if this_cutoff is 0 or NA (indicating a pseudo taxa)
  if (is.na(this_cutoff) || this_cutoff == 0) {
    # Get cutoff for pseudo taxa from pseudo_cutoffs
    this_cutoff <- pseudo_cutoffs %>%
      filter(rank == this_rank) %>%
      pull(cutoff)
    
    message(paste0("Assigning cutoff ", this_cutoff, " to pseudo taxa ", this_supertaxon))
  }
  
  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
    
    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }
    
    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }
    
    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]
    
    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")
    
    ## Cluster unidentified ASVs ####
    
    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))
    
    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )
    
    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length bitscore qstart qend'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )
    
    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        coverage = V4,
        bitscore = V5,
        start = V6,
        end = V7
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score with the necessary adjustment
      mutate(
        score = sim / 100,
        score = case_when(
          coverage < minlen ~ (score * coverage) / minlen,
          TRUE ~ score
        )
      ) %>%
      filter(score > this_cutoff) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        class = this_supertaxon,
        family = "unidentified",
        genus = "unidentified",
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        confidence = 0,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        class, order, family, genus, species,
        rank, cutoff, confidence, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      )
    
    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )
    
    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)
    
    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
    
    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }
  
  ### (4b) De novo clustering ###########################################
  
  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)
  
  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))
    
    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")
    
    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))
    
    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")
    
    # Initialise pseudo_id
    pseudo_id <- 1
    
    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)
    
    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (5) Family clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "order"
this_rank <- "family"
this_subrank <- "genus"

# Write the order file as the family file
fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Order cut-offs
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (5a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]
  
  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))
  
  # Check if this_cutoff is 0 or NA (indicating a pseudo taxa)
  if (is.na(this_cutoff) || this_cutoff == 0) {
    # Get cutoff for pseudo taxa from pseudo_cutoffs
    this_cutoff <- pseudo_cutoffs %>%
      filter(rank == this_rank) %>%
      pull(cutoff)
    
    message(paste0("Assigning cutoff ", this_cutoff, " to pseudo taxa ", this_supertaxon))
  }
  
  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
    
    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }
    
    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }
    
    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]
    
    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")
    
    ## Cluster unidentified ASVs ####
    
    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))
    
    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )
    
    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length bitscore qstart qend'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )
    
    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        coverage = V4,
        bitscore = V5,
        start = V6,
        end = V7
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score with the necessary adjustment
      mutate(
        score = sim / 100,
        score = case_when(
          coverage < minlen ~ (score * coverage) / minlen,
          TRUE ~ score
        )
      ) %>%
      filter(score > this_cutoff) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        order = this_supertaxon,
        genus = "unidentified",
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        confidence = 0,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        order, family, genus, species,
        rank, cutoff, confidence, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      )
    
    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )
    
    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)
    
    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
    
    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }
  
  ### (5b) De novo clustering ###########################################
  
  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)
  
  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))
    
    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")
    
    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))
    
    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")
    
    # Initialise pseudo_id
    pseudo_id <- 1
    
    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)
    
    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (6) Genus clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "family"
this_rank <- "genus"
this_subrank <- "species"

# Write the family file as the genus file
fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Family cut-offs
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (6a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]
  
  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))
  
  # Check if this_cutoff is 0 or NA (indicating a pseudo taxa)
  if (is.na(this_cutoff) || this_cutoff == 0) {
    # Get cutoff for pseudo taxa from pseudo_cutoffs
    this_cutoff <- pseudo_cutoffs %>%
      filter(rank == this_rank) %>%
      pull(cutoff)
    
    message(paste0("Assigning cutoff ", this_cutoff, " to pseudo taxa ", this_supertaxon))
  }
  
  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
    
    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }
    
    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }
    
    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]
    
    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")
    
    ## Cluster unidentified ASVs ####
    
    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))
    
    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )
    
    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length bitscore qstart qend'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )
    
    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        coverage = V4,
        bitscore = V5,
        start = V6,
        end = V7
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score with the necessary adjustment
      mutate(
        score = sim / 100,
        score = case_when(
          coverage < minlen ~ (score * coverage) / minlen,
          TRUE ~ score
        )
      ) %>%
      filter(score > this_cutoff) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        family = this_supertaxon,
        species = "unidentified",
        rank = this_rank,
        cutoff = this_cutoff,
        confidence = 0,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff"),
          paste0(this_subrank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        family, genus, species,
        rank, cutoff, confidence, score,
        paste0(this_rank, "_cutoff"),
        paste0(this_subrank, "_cutoff")
      )
    
    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )
    
    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)
    
    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
    
    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }
  
  ### (6b) De novo clustering ###########################################
  
  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)
  
  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))
    
    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")
    
    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))
    
    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")
    
    # Initialise pseudo_id
    pseudo_id <- 1
    
    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)
    
    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

###############################################################################*
# (7) Species clusters ##########################################################
###############################################################################*

## Define the ranks of interest and its enclosing supertaxon
this_superrank <- "genus"
this_rank <- "species"

# Write the genus file as the species file
fread(paste0("./tmp_clusters/", this_superrank, "_clusters.txt")) %>%
  fwrite(., paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))

# Genus cut-offs
supertaxa_cutoffs <- taxa_cutoffs %>%
  select(paste(this_superrank), paste0(this_rank, "_cutoff")) %>%
  filter(
    !get(this_superrank) %in% c("unidentified", "")
  ) %>%
  unique(.)

message(paste0("Starting clustering for ", this_rank, " !!!"))

### (7a) Reference-based clustering ############################################

# Loop over each supertaxon and apply the respective cutoff
for (i in 1:nrow(supertaxa_cutoffs)) {
  taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
  this_supertaxon <- supertaxa_cutoffs[[this_superrank]][i]
  this_cutoff <- supertaxa_cutoffs[[paste0(this_rank, "_cutoff")]][i]
  
  message(paste0("Starting clustering for ", this_supertaxon, " at rank ", this_rank, " using cutoff ", this_cutoff, "..."))
  
  # Check if this_cutoff is 0 or NA (indicating a pseudo taxa)
  if (is.na(this_cutoff) || this_cutoff == 0) {
    # Get cutoff for pseudo taxa from pseudo_cutoffs
    this_cutoff <- pseudo_cutoffs %>%
      filter(rank == this_rank) %>%
      pull(cutoff)
    
    message(paste0("Assigning cutoff ", this_cutoff, " to pseudo taxa ", this_supertaxon))
  }
  
  repeat {
    #### Create cluster cores ####
    taxa_cutoffs <- fread(paste0("./tmp_clusters/", this_rank, "_clusters.txt"))
    
    identified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) != "unidentified" & get(this_rank) != "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no identified ASVs
    initial_identified_count <- length(identified_asvs)
    if (initial_identified_count == 0) {
      break
    }
    
    unidentified_asvs <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.)
    
    # Skip the loop for a given taxon if there are no unidentified ASVs
    initial_unidentified_count <- length(unidentified_asvs)
    if (initial_unidentified_count == 0) {
      break
    }
    
    # Filter identified and unidentified sequences
    identified_sequences <- asv_sequences[names(asv_sequences) %in% identified_asvs]
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% unidentified_asvs]
    
    # Write the fasta files
    writeXStringSet(identified_sequences, "./tmp/identified_fasta")
    writeXStringSet(unidentified_sequences, "./tmp/unidentified_fasta")
    
    ## Cluster unidentified ASVs ####
    
    message(paste0("Clustering unidentified ASVs to ", this_supertaxon, "..."))
    
    # Build the BLAST database using identified ASVs
    system2("makeblastdb",
            args = c(
              "-in", "./tmp/identified_fasta",
              "-dbtype", "nucl"
            )
    )
    
    # BLAST unidentified ASVs against the cluster cores
    system2("blastn",
            args = c(
              "-query", "./tmp/unidentified_fasta",
              "-db", "./tmp/identified_fasta",
              "-outfmt", "'6 qseqid sseqid pident length bitscore qstart qend'",
              "-task", "blastn-short",
              "-num_threads", as.character(threads),
              "-out", "./tmp/unidentified_blast.out"
            )
    )
    
    # Read and filter BLAST results, and update taxonomy
    new_clusters <- fread("./tmp/unidentified_blast.out") %>%
      select(
        OTU_ID = V1,
        reference_ID = V2,
        sim = V3,
        coverage = V4,
        bitscore = V5,
        start = V6,
        end = V7
      ) %>%
      group_by(OTU_ID) %>%
      dplyr::slice(1) %>%
      ungroup(.) %>%
      # Compute score with the necessary adjustment
      mutate(
        score = sim / 100,
        score = case_when(
          coverage < minlen ~ (score * coverage) / minlen,
          TRUE ~ score
        )
      ) %>%
      filter(score > this_cutoff) %>%
      # Initialise taxonomy columns with placeholders (including dynamic rank)
      mutate(
        genus = this_supertaxon,
        rank = this_rank,
        cutoff = this_cutoff,
        confidence = 0,
        score = score
      ) %>%
      # Join based on reference_ID in new assignments
      left_join(
        taxa_cutoffs %>% select(
          reference_ID = OTU_ID, !!sym(this_rank),
          paste0(this_rank, "_cutoff")
        ),
        by = "reference_ID"
      ) %>%
      # Select only the relevant columns for the output
      select(
        OTU_ID, reference_ID,
        genus, species,
        rank, cutoff, confidence, score,
        paste0(this_rank, "_cutoff")
      )
    
    # Update taxonomy
    taxa_cutoffs <- bind_rows(
      new_clusters,
      taxa_cutoffs %>% filter(!OTU_ID %in% new_clusters[["OTU_ID"]])
    )
    
    # Check if the number of unidentified ASVs has changed
    remaining_unidentified_count <- taxa_cutoffs %>%
      filter(get(this_superrank) == this_supertaxon) %>%
      filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
      select(OTU_ID) %>%
      pull(.) %>%
      length(.)
    
    # Write the updated taxonomy
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
    
    ## Break or repeat the process #####
    # Break the loop if there are no more unidentified ASVs
    if (remaining_unidentified_count == 0) {
      message(paste0("There are no more unidentified ASVs for ", this_supertaxon, "..."))
      break
      # Break the loop if no new assignments on a repeated clustering round
    } else if (remaining_unidentified_count == initial_unidentified_count) {
      message(paste0("No new assignments were made for ", this_supertaxon, "..."))
      break
      # Otherwise, continue clustering until one of the above conditions is met
    } else {
      message(paste0("Continuing clustering with remaining unidentified ASVs in ", this_supertaxon, "..."))
    }
  }
  
  ### (7b) De novo clustering ###########################################
  
  # Remove ASVs that were already clustered in the reference-based clustering
  remaining_unidentified_asvs <- taxa_cutoffs %>%
    filter(get(this_superrank) == this_supertaxon) %>%
    filter(get(this_rank) == "unidentified" | get(this_rank) == "") %>%
    select(OTU_ID) %>%
    pull(.)
  
  # If no remaining unidentified ASVs, skip the de novo clustering step
  if (length(remaining_unidentified_asvs) == 0) {
    message(paste0("No remaining ASVs for de novo clustering of ", this_supertaxon, "..."))
  } else {
    message(paste0("Commencing de novo clustering of unidentified ASVs in ", this_supertaxon, "..."))
    
    # Write unidentified sequences to FASTA for de novo clustering
    unidentified_sequences <- asv_sequences[names(asv_sequences) %in% remaining_unidentified_asvs]
    writeXStringSet(unidentified_sequences, "./tmp/remaining_unidentified.fasta")
    
    # Perform de novo clustering using blastclust
    identity_cutoff <- this_cutoff * 100
    system2("blastclust", args = c(
      "-i", paste0("./tmp/remaining_unidentified.fasta"),
      "-S", as.character(identity_cutoff),
      "-a", as.character(threads),
      "-p", "F",
      "-o", paste0("./tmp/de_novo_clusters.txt")
    ))
    
    # Read the de novo clustering output as a single column
    pseudo_clusters <- fread("./tmp/de_novo_clusters.txt", header = FALSE, sep = "\n", col.names = "cluster")
    
    # Initialise pseudo_id
    pseudo_id <- 1
    
    # Process the clusters and update the taxonomy file
    taxa_cutoffs <- pseudo_clusters %>%
      # Split each row (cluster) into individual ASVs
      mutate(ASVs = str_split(cluster, " ")) %>%
      # Expand the ASVs (turn each list into individual rows)
      unnest(ASVs) %>%
      # Group by the original cluster
      group_by(cluster) %>%
      # Assign a unique pseudo cluster name for each cluster using cur_group_id()
      mutate(
        pseudo_name = paste0(this_supertaxon, "_pseudo_", this_rank, "_", sprintf("%04d", cur_group_id()))
      ) %>%
      ungroup() %>%
      # Retain the OTU ID with the size information
      mutate(OTU_ID = ASVs) %>%
      select(OTU_ID, pseudo_name) %>%
      # Join with existing taxa_cutoffs on OTU_ID
      full_join(
        taxa_cutoffs,
        by = "OTU_ID"
      ) %>%
      # Update the this_rank column with the pseudo_name for the new clusters
      mutate(
        !!sym(this_rank) := case_when(
          !is.na(pseudo_name) ~ pseudo_name,
          TRUE ~ !!sym(this_rank)
        )
      ) %>%
      # Remove the pseudo_name column after updating
      select(-pseudo_name)
    
    # Write the updated taxonomy file
    fwrite(taxa_cutoffs, paste0("./tmp_clusters/", this_rank, "_clusters.txt"), sep = "\t")
  }
}

message(paste0("Clustering for ", this_rank, " complete!!!"))

# Empty the temporary folder
unlink(files, recursive = TRUE)

# Clean up the environment
suppressWarnings({
  rm(list = c(
    ls(pattern = "this_"),
    ls(pattern = "unidentified_"),
    ls(pattern = "identified_"),
    ls(pattern = "new_")
  ))
})

EOF