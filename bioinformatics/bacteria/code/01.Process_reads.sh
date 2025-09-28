#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Define directories
export INPUT="../data/raw_data"
export OUTPUT="../data/"
export PRIMERS_TRIMMED_DIR="$OUTPUT/tmp_primers_trimmed"

# Make directories
mkdir -p "$OUTPUT" "$PRIMERS_TRIMMED_DIR"

# Define constants
readonly PRIMER_FWD="GTGYCAGCMGCCGCGGTAA"    # Forward primer 515F (Parada) for the V4 region
readonly PRIMER_REV="GGACTACNVGGGTWTCTAAT"   # Reverse primers 806R (Apprill) for the V4 region
readonly MIN_OVERLAP_FWD=17
readonly MIN_OVERLAP_REV=20
readonly MAX_ERROR_RATE=2

# Activate the conda environment
source ~/.bashrc
conda activate dynamic_clustering

################################################################################
# Prepare reads for DADA2
################################################################################

trim_primers() {
    mkdir -p "$PRIMERS_TRIMMED_DIR" "$OUTPUT"

    echo "Trimming primers from forward reads..."
    for file in "${INPUT}"/*_R1_001.fastq.gz; do
        base_name=$(basename "$file" | cut -d'_' -f1)
        cutadapt \
            -g "${PRIMER_FWD};min_overlap=${MIN_OVERLAP_FWD}" \
            -o "${PRIMERS_TRIMMED_DIR}/${base_name}_R1.fastq.gz" "$file" \
            -e "$MAX_ERROR_RATE" \
            --times 4 \
            --discard-untrimmed >> "${OUTPUT}/logfile_cutadapt.txt" 2>&1
    done

    echo "Trimming primers from reverse reads..."
    for file in "${INPUT}"/*_R2_001.fastq.gz; do
        base_name=$(basename "$file" | cut -d'_' -f1)
        cutadapt \
            -g "${PRIMER_REV};min_overlap=${MIN_OVERLAP_REV}" \
            -o "${PRIMERS_TRIMMED_DIR}/${base_name}_R2.fastq.gz" "$file" \
            -e "$MAX_ERROR_RATE" \
            --times 4 \
            --discard-untrimmed >> "${OUTPUT}/logfile_cutadapt.txt" 2>&1
    done

    echo "Synchronizing reads after trimming..."
    for fwd_file in "${PRIMERS_TRIMMED_DIR}"/*_R1.fastq.gz; do
        rev_file="${fwd_file/_R1/_R2}"
        seqkit pair -1 "$fwd_file" -2 "$rev_file" -O "${OUTPUT}" >> "${OUTPUT}/logfile_seqkit.txt" 2>&1
    done
}

# Execute the primer trimming function
trim_primers

################################################################################
# Run DADA2
################################################################################

Rscript --vanilla - <<EOF
library(dada2); packageVersion("dada2")
library(Biostrings)
library(tidyverse)

path <- "../data" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), \`[\`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(175,125),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

# Learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Track reads through the pipeline
#getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#rownames(track) <- sample.names
#data.table::fwrite(track, file.path(path, "track.txt"), sep="\t")

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, file.path(path, "silva_nr99_v138.1_wSpecies_train_set.fa.gz"), multithread=TRUE)

# Save taxonomy and ASV tables
seqtab.nochim %>% 
  t() %>% 
  as.data.frame(stringsAsFactors = F) %>% 
  rownames_to_column("ASV_ID") %>%
  data.table::fwrite(., file.path(path, "ASVs.txt"), sep="\t")
taxa %>%
  as.data.frame(stringsAsFactors = T) %>%
  rownames_to_column("ASV_ID") %>%
  data.table::fwrite(., file.path(path, "taxonomy.txt"), sep="\t")

# Save ASV sequences as FASTA
asv_seqs <- colnames(seqtab.nochim)
asv_ids <- paste0("OTU_", seq_len(length(asv_seqs)))
names(asv_seqs) <- asv_ids
asv_stringset <- DNAStringSet(asv_seqs)
writeXStringSet(asv_stringset, file.path(path, "ASVs.fasta"))

EOF

# Deactivate the conda environemtn
conda deactivate