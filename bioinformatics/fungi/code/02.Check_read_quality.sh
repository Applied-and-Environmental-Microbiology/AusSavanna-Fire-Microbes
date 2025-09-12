#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=7-00:00:00
#SBATCH --partition=month
#SBATCH --output=slurm/%x.%j.out

# Script: Remove primers, extract the ITS region of interest, and quality filter reads

# Software:
# --------
# ITSxpress v2.1.0: https://github.com/USDA-ARS-GBRU/itsxpress
# Cutadapt v4.9: https://cutadapt.readthedocs.io/
# VSEARCH v2.22.1: https://github.com/torognes/vsearch
# SeqKit v2.8.2: https://github.com/rcedgar/usearch12
# FastQC v0.12.1: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# MultiQC v1.24.1: https://multiqc.info/

# Note: Combined lentgth of 370 for forward and reverse reads will cover > 99%
# of unique ITS sequences in UNITE. The minimum overlap for DADA2 is 12.


# Define directories
export RAW_DATA_DIR="../data/raw_data"

# Activate the conda environment
echo "Starting quality check at: $(date)"
source ~/.bashrc
conda activate sequence_prep

###################################################
# Quality Check
###################################################

echo "Checking quality of raw data..."

# Run FastQC
fastqc \
    -t "$NUM_THREADS" \
    -o "$RAW_DATA_DIR" \
    "${RAW_DATA_DIR}"/*_R1_001.fastq.gz \
    "${RAW_DATA_DIR}"/*_R2_001.fastq.gz

# Combine FastQC reports into a single report
multiqc "$RAW_DATA_DIR"/. \
    -o "$RAW_DATA_DIR"

# Clean up FastQC output - delete individual reports
rm "${RAW_DATA_DIR}"/*fastqc.html "${RAW_DATA_DIR}"/*fastqc.zip
rm -r "${RAW_DATA_DIR}"/multiqc_report_data
mv "${RAW_DATA_DIR}"/multiqc_report.html "${RAW_DATA_DIR}"/quality_report.html

conda deactivate
echo "Finished quality check at: $(date)"