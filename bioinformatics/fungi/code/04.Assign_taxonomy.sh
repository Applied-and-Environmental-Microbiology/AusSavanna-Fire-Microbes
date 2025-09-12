#!/bin/bash -e
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

## Script: Taxonomic assignment with dnabarcoder
## Credit: https://github.com/vuthuyduong/dnabarcoder
## Date: 16th July 2024
## Software:
##  - BLAST v2.14.1 - https://blast.ncbi.nlm.nih.gov/Blast.cgi
##  - dnabarcoder v1.0.6 - https://github.com/vuthuyduong/dnabarcoder

# Constants and subdirectories
readonly REFERENCE_UNIQUES="dnabarcoder/unite2024ITS2.unique.fasta"
readonly CLASSIFIER="dnabarcoder/unite2024ITS2.classification"
readonly BEST_MATCH="dnabarcoder/unite2024ITS2.unique.cutoffs.best.json"
readonly QUERY_SEQUENCES="../data/prepared_reads/ASVs.fasta"
readonly OUTPUT="../data/taxonomy/"

declare -A REFERENCE_UNIQUE_RANKS=(
    ["species"]="dnabarcoder/unite2024ITS2.unique.species.fasta"
    ["genus"]="dnabarcoder/unite2024ITS2.unique.genus.fasta"
    ["family"]="dnabarcoder/unite2024ITS2.unique.family.fasta"
    ["order"]="dnabarcoder/unite2024ITS2.unique.order.fasta"
    ["class"]="dnabarcoder/unite2024ITS2.unique.class.fasta"
    ["phylum"]="dnabarcoder/unite2024ITS2.unique.phylum.fasta"
)

################################################################################
# Main script ##################################################################
################################################################################

echo "Starting at: $(date)"

# Activate the conda environment
source ~/.bashrc
conda activate dynamic_clustering

# Reduce the size of the reference dataset by retaining unique sequences
if ! python dnabarcoder/aidscripts/selectsequences.py \
    -i "$REFERENCE_SEQUENCES" \
    -unique yes \
    -o "$REFERENCE_UNIQUES"; then
    echo "Error: Failed to reduce reference dataset" >&2
    exit 1
fi

# Process each taxonomic rank
for rank in "${!REFERENCE_UNIQUE_RANKS[@]}"; do
    output="${REFERENCE_UNIQUE_RANKS[$rank]}"
    
    # Select sequences by rank
    if ! python dnabarcoder/aidscripts/selectsequences.py \
        -i "$REFERENCE_UNIQUES" \
        -c "$CLASSIFIER" \
        -rank "$rank" \
        -o "$output"; then
        echo "Error: Failed to select sequences for rank: $rank" >&2
        exit 1
    fi
    
    # Search for the best matches of the sequences
    if ! python dnabarcoder/dnabarcoder.py search \
        -i "$QUERY_SEQUENCES" \
        -r "$output" \
        -ml 50; then
        echo "Error: Failed to search for best matches for rank: $rank" >&2
        exit 1
    fi
    
    # Classify the sequences to different taxonomic groups
    if ! python dnabarcoder/dnabarcoder.py classify \
        -i "dnabarcoder/ASVs.unite2024ITS2.unique.${rank}_BLAST.bestmatch" \
        -c "$CLASSIFIER" \
        -cutoffs "$BEST_MATCH"; then
        echo "Error: Failed to classify sequences for rank: $rank" >&2
        exit 1
    fi
    
    # Move results to output directory
    mv "dnabarcoder/ASVs.unite2024ITS2.unique.${rank}_BLAST.classification" "$OUTPUT/${rank}_classification.txt"
done

# Merge classification files giving prefernce to lower ranks
Rscript -e "taxa_dir <- '${OUTPUT}'; source('functions/merge_taxonomy.R')"

# Deactivate conda environment
conda deactivate

echo "Finished at: $(date)"
