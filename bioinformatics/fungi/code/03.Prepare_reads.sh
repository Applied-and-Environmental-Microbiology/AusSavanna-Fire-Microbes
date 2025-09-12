#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=slurm/%x.%j.out

# Script: Remove primers, quality filter, and denoise reads

# Software:
# --------
# R v3.3.1 -- Beagle Scouts: https://www.r-project.org/
# tidyverse v2.0.0: https://tidyverse.tidyverse.org
# DADA2 v1.30.0: https://benjjneb.github.io/dada2/
# Cutadapt v4.9: https://cutadapt.readthedocs.io/
# VSEARCH v2.22.1: https://github.com/torognes/vsearch
# SeqKit v2.8.2: https://bioinf.shenwei.me/seqkit/
# FastQC v0.12.1: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# MultiQC v1.24.1: https://multiqc.info/

# Note: Combined lentgth of 370 for forward and reverse reads will cover > 99%
# of unique ITS sequences in UNITE. The minimum overlap for DADA2 is 12:
# 221 + 161 - 370 = 12

# Source the functions file
source prepare_reads_functions.sh

# Define directories
export INPUT="../data/raw_data"
export OUTPUT="../data/prepared_reads"

# Define constants
export NUM_THREADS=40
export PRIMER_FWD="GTGARTCATCGAATCTTTG" # Cutadapt
export PRIMER_REV="TCCTCCGCTTATTGATATGC" # Cutadapt
export PRIMER_FWD_RC=$(RC "$PRIMER_FWD") # Cutadapt
export PRIMER_REV_RC=$(RC "$PRIMER_REV") # Cutadapt
export MIN_OVERLAP_FWD=16 # Cutadapt
export MIN_OVERLAP_REV=17 # Cutadapt
export MAX_ERROR_RATE=3 # Cutadapt
export MAXEE=2 # VSEARCH
export MINLEN=50 # VSEARCH
export MAXLEN_FWD=221 # VSEARCH
export MAXLEN_REV=161 # VSEARCH
export MAXN=0 # VSEARCH
export QMAX=41 # VSEARCH
export PRECLUSTERED_IDENTITY=0.98 # VSEARCH
export REFERENCE_SEQS="../data/reference_data/unite2024ITS2.all.fasta" # VSEARCH
export REGION="ITS2"  # ITSx: The target ITS region: SSU, ITS1, 5.8S, ITS2, LSU, all, none. Outputs only the actual ITS sequences (ITS1, ITS2) by default.
export TAXA="Fungi"   # ITSx: The taxonomic group to extract: ALL, Fungi, Metazoa, Viridiplantae, Streptophyta, Rhodophyta, Protozoa, Chromista, Bacteria, Archaea, none. Extracts all taxa by default.

# Make directories
mkdir -p "$OUTPUT"

# Activate the conda environment
echo "Starting at:" $(date)
source ~/.bashrc
#conda activate sequence_prep
conda activate dynamic_clustering

trim_primers
qual_filter
qual_check
denoise
chimera_filter
extract_ITS
generate_ASV_table

# Cleanup
#rm -r "$PRIMERS_TRIMMED_DIR" "$SYNC_PAIRS_DIR" "$QUAL_FILTERED_DIR" "$SYNC_FILTERED_DIR" "$DENOISED_DIR" "$CHIMERA_FILTERED_DIR"

# Deactive the conda environment
conda deactivate
echo "Finished at:" $(date)
