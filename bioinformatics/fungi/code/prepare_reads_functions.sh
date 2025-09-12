
export PRIMERS_TRIMMED_DIR="${OUTPUT}/tmp_primers_trimmed"
export SYNC_PAIRS_DIR="${OUTPUT}/tmp_sync_primers_trimmed"
export QUAL_FILTERED_DIR="${OUTPUT}/tmp_qual_filtered"
export SYNC_FILTERED_DIR="${OUTPUT}/tmp_sync_qual_filtered"
export DENOISED_DIR="${OUTPUT}/tmp_denoised"
export CHIMERA_FILTERED_DIR="${OUTPUT}/tmp_chimera_filtered"
export ITSX_DIR="${OUTPUT}/tmp_itsx"

################################################################################
# Function to construct reverse-complement sequences
################################################################################
RC() {
    echo "$1" | tr "ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv" "TACGAtacgaNnRrYySsWwMmKkVvHhDdBb" | rev
}

################################################################################
# Function to trim primers using Cutadapt
################################################################################
trim_primers() {
    mkdir -p "$PRIMERS_TRIMMED_DIR" "$SYNC_PAIRS_DIR"

    echo "Trimming primers from forward reads..."
    for file in "${INPUT}"/*_R1_001.fastq.gz; do
        base_name=$(basename "$file" | cut -d'_' -f1)
        cutadapt \
            -g "${PRIMER_FWD};required;min_overlap=${MIN_OVERLAP_FWD}"..."${PRIMER_REV_RC};optional;min_overlap=${MIN_OVERLAP_REV}" \
            -o "${PRIMERS_TRIMMED_DIR}/${base_name}_R1.fastq.gz" "$file" \
            -e "$MAX_ERROR_RATE" \
            --times 2 \
            --discard-untrimmed >> "${OUTPUT}/logfile_cutadapt.txt" 2>&1
    done

    echo "Trimming primers from reverse reads..."
    for file in "${INPUT}"/*_R2_001.fastq.gz; do
        base_name=$(basename "$file" | cut -d'_' -f1)
        cutadapt \
            -g "${PRIMER_REV};required;min_overlap=${MIN_OVERLAP_REV}"..."${PRIMER_FWD_RC};optional;min_overlap=${MIN_OVERLAP_FWD}" \
            -o "${PRIMERS_TRIMMED_DIR}/${base_name}_R2.fastq.gz" "$file" \
            -e "$MAX_ERROR_RATE" \
            --times 2 \
            --discard-untrimmed >> "${OUTPUT}/logfile_cutadapt.txt" 2>&1
    done

    echo "Synchronizing reads after trimming..."
    for fwd_file in "${PRIMERS_TRIMMED_DIR}"/*_R1.fastq.gz; do
        rev_file="${fwd_file/_R1/_R2}"
        seqkit pair -1 "$fwd_file" -2 "$rev_file" -O "${SYNC_PAIRS_DIR}" >> "${OUTPUT}/logfile_seqkit_1.txt" 2>&1
    done
}

################################################################################
# Function to quality filter reads using VSEARCH
################################################################################
qual_filter() {
    mkdir -p "$QUAL_FILTERED_DIR" "$SYNC_FILTERED_DIR"

    echo "Quality filtering forward reads..."
    for file in "${SYNC_PAIRS_DIR}"/*_R1.fastq.gz; do
        vsearch \
            --fastq_filter "$file" \
            --fastq_maxee "$MAXEE" \
            --fastq_minlen "$MINLEN" \
            --fastq_trunclen_keep "$MAXLEN_FWD" \
            --fastq_maxns "$MAXN" \
            --fastq_qmax "$QMAX" \
            --fastqout "${QUAL_FILTERED_DIR}/$(basename "$file" .fastq.gz).fastq" 2>> "${OUTPUT}/logfile_vsearch.txt"
    done

    echo "Quality filtering reverse reads..."
    for file in "${SYNC_PAIRS_DIR}"/*_R2.fastq.gz; do
        vsearch \
            --fastq_filter "$file" \
            --fastq_maxee "$MAXEE" \
            --fastq_minlen "$MINLEN" \
            --fastq_trunclen_keep "$MAXLEN_REV" \
            --fastq_maxns "$MAXN" \
            --fastq_qmax "$QMAX" \
            --fastqout "${QUAL_FILTERED_DIR}/$(basename "$file" .fastq.gz).fastq" 2>> "${OUTPUT}/logfile_vsearch.txt"
    done
    pigz -p "$NUM_THREADS" "${QUAL_FILTERED_DIR}"/*.fastq

    echo "Synchronizing reads after quality filtering..."
    for fwd_file in "${QUAL_FILTERED_DIR}"/*_R1.fastq.gz; do
        rev_file="${fwd_file/_R1/_R2}"
        seqkit pair -1 "$fwd_file" -2 "$rev_file" -O "${SYNC_FILTERED_DIR}" >> "${OUTPUT}/logfile_seqkit_2.txt" 2>&1
    done
}

################################################################################
# Function for quality check using FastQC and MultiQC
################################################################################
qual_check() {
    echo "Checking quality of quality filtered reads..."
    fastqc -t "$NUM_THREADS" -o "$OUTPUT" "${SYNC_FILTERED_DIR}"/*fastq.gz
    multiqc "$OUTPUT" -o "$OUTPUT"
    rm "${OUTPUT}"/*fastqc.html "${OUTPUT}"/*fastqc.zip
    rm -r "${OUTPUT}/multiqc_data"
    mv "${OUTPUT}/multiqc_report.html" "${OUTPUT}/quality_filtered_report.html"
}

################################################################################
# Function to denoise using DADA2 in R
################################################################################
denoise() {
    mkdir -p "$DENOISED_DIR" "$CHIMERA_FILTERED_DIR"
    echo "Denoising reads using DADA2..."
    Rscript - <<EOF
        suppressPackageStartupMessages(require(dada2))
        suppressPackageStartupMessages(require(tidyverse))
        filtFs <- sort(list.files("$SYNC_FILTERED_DIR", pattern="R1.fastq.gz", full.names = TRUE))
        filtRs <- sort(list.files("$SYNC_FILTERED_DIR", pattern="R2.fastq.gz", full.names = TRUE))
        sample.names <- sapply(strsplit(basename(filtFs), "_"), \`[\`, 1)
        names(filtFs) <- sample.names
        names(filtRs) <- sample.names
        errF <- learnErrors(filtFs, multithread=8)
        errR <- learnErrors(filtRs, multithread=8)
        dadaFs <- dada(filtFs, err=errF, multithread=8)
        dadaRs <- dada(filtRs, err=errR, multithread=8)
        mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
        seqtab <- makeSequenceTable(mergers)
        getN <- function(x) sum(getUniques(x))
        track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab))
        saveRDS(seqtab, file.path("$DENOISED_DIR", "seqtab.rds"))
        write.table(track, file.path("$OUTPUT", "logfile_dada2.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

        # With chimera removal
        #seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=8, verbose=TRUE)
        #track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
        #saveRDS(seqtab.nochim, file.path("$DENOISED_DIR", "seqtab.rds"))

        # .rds to .txt with VSEARCH formatted headers
        as.data.frame(seqtab) %>%
            rownames_to_column('sample.names') %>%
            pivot_longer(-sample.names, names_to = 'seq', values_to = 'size') %>%
            filter(size > 0) %>%
            ungroup() %>%
            mutate(seq.name = paste0(sample.names, ".", row_number(), ";size=", size)) %>%
            select(seq.name, seq) %>%
            write.table(., file = file.path("$DENOISED_DIR", "seqtab.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
EOF
    if [ $? -ne 0 ]; then
        echo "R script failed"
        exit 1
    fi

    seqkit tab2fx $DENOISED_DIR/seqtab.txt -o $CHIMERA_FILTERED_DIR/all.fasta
}

## Function for dereplication across samples, de novo and referenced-based chimera filtering
chimera_filter() {
    local files_to_remove=("all.denovo.nonchimeras.fasta" "all.derep.fasta" "all.derep.uc" "all.nonchimeras.derep.fasta" "all.nonchimeras.fasta" "all.preclustered.fasta" "all.preclustered.uc" "all.ref.nonchimeras.fasta")

    for file in "${files_to_remove[@]}"; do
        rm -f "$CHIMERA_FILTERED_DIR/$file"
    done

    echo 'Dereplicating across samples...'
    vsearch \
        --derep_fulllength "$CHIMERA_FILTERED_DIR/all.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.derep.uc" \
        --output "$CHIMERA_FILTERED_DIR/all.derep.fasta"

    echo 'Preclustering reads...'
    vsearch \
        --cluster_size "$CHIMERA_FILTERED_DIR/all.derep.fasta" \
        --threads $NUM_THREADS \
        --id $PRECLUSTERED_IDENTITY \
        --strand plus \
        --sizein --sizeout \
        --fasta_width 0 \
        --uc "$CHIMERA_FILTERED_DIR/all.preclustered.uc" \
        --centroids "$CHIMERA_FILTERED_DIR/all.preclustered.fasta"

    echo 'Starting de novo chimera detection...'
    vsearch \
        --uchime3_denovo "$CHIMERA_FILTERED_DIR/all.preclustered.fasta" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta"

    echo 'Starting reference-based chimera detection...'
    vsearch \
        --uchime_ref "$CHIMERA_FILTERED_DIR/all.denovo.nonchimeras.fasta" \
        --threads "$NUM_THREADS" \
        --db "$REFERENCE_SEQS" \
        --sizein --sizeout \
        --fasta_width 0 \
        --nonchimeras "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta"

    echo 'Extracting all non-chimeric sequences...'
    perl map.pl "$CHIMERA_FILTERED_DIR/all.derep.fasta" "$CHIMERA_FILTERED_DIR/all.preclustered.uc" "$CHIMERA_FILTERED_DIR/all.ref.nonchimeras.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta"
    #perl map.pl "$CHIMERA_FILTERED_DIR/all.fasta" "$CHIMERA_FILTERED_DIR/all.derep.uc" "$CHIMERA_FILTERED_DIR/all.nonchimeras.derep.fasta" > "$CHIMERA_FILTERED_DIR/all.nonchimeras.fasta"
}

extract_ITS() { 
    mkdir -p "$ITSX_DIR"
    if [ ! -s "${CHIMERA_FILTERED_DIR}/all.nonchimeras.derep.fasta" ]; then
        echo "Error: Input file all.nonchimeras.derep.fasta is empty!"
        exit 1
    fi

    ITSx -i "${CHIMERA_FILTERED_DIR}/all.nonchimeras.derep.fasta" \
        --save_regions $REGION -t F \
        —cpu "$NUM_THREADS" -o "${ITSX_DIR}/ITSX"
}

## Function for generating an ASV table
generate_ASV_table() {

    # Rename ITSx headers
    seqkit replace -p '\|.*' -r '' "${ITSX_DIR}/ITSX.ITS2.fasta" -o "${ITSX_DIR}/VSEARCH.ITS2.fasta"

    # Map derep info
    perl map.pl "${CHIMERA_FILTERED_DIR}/all.fasta" "${CHIMERA_FILTERED_DIR}/all.derep.uc" "${ITSX_DIR}/VSEARCH.ITS2.fasta" > "${OUTPUT}/all.ITS2.fasta"
    
    # Generate ASV table
    echo "Generating ASV table at..."

    vsearch \
        --cluster_unoise "${OUTPUT}/all.ITS2.fasta" \
        --threads $NUM_THREADS \
        --sizein --sizeout \
        --relabel ASV_ \
        --fasta_width 0 \
        --uc "${OUTPUT}/ASVs.uc" \
        --biomout "${OUTPUT}/ASVs.biom" \
        --otutabout "${OUTPUT}/ASVs.txt" \
        --centroids "${OUTPUT}/ASVs.fasta" 2>> "${OUTPUT}/logfile_ASVs.txt"

    # Rename the header in the file
    sed -i '1s/#OTU ID/OTU_ID/' "${OUTPUT}/ASVs.txt"
}