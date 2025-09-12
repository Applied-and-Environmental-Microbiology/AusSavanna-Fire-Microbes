
# Activate the environment
conda activate gemelli_env

# Create a directory for distance data
mkdir data/distances

# Convert .txt to .biom
biom convert -i output/otus_prevalence_filtered.txt -o data/distances/otus.biom --table-type="OTU table" --to-json

# Run gemelli based on identity
gemelli rpca --in-biom data/distances/otus.biom --output-dir data/distances

# Check the output names
ls data/distances

# Move the distance matrix to output directory
mv data/distances/distance-matrix.tsv output/distance_fungi.txt
conda run -n gemelli_env biom convert -i output/otus_prevalence_filtered.txt -o data/distances/otus.biom --table-type="OTU table" --to-json