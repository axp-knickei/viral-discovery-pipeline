#!/bin/bash

# Viro-Flow Stage 1: Per-Sample Viral Discovery Pipeline
#
# This script processes a single sequencing sample (paired-end or single-end)
# from raw FASTQ reads to high-quality, host-removed viral contigs.
# For more details, see the README.md file.

# --- Script Setup ---
# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error when substituting.
set -u
# Pipes will fail if any command in the pipeline fails.
set -o pipefail

# --- Logging Function ---
# Usage: log "Your message here"
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] -" "$1"
}

# --- Usage Function ---
usage() {
    echo "Usage: $0 -s <SAMPLE_ID> [-1 <READ1_FQ> -2 <READ2_FQ> | -r <SINGLE_FQ>]"
    echo "  -s  SAMPLE_ID      A unique identifier for the sample."
    echo "  -1  READ1_FQ       Path to the forward/R1 FASTQ file (paired-end)."
    echo "  -2  READ2_FQ       Path to the reverse/R2 FASTQ file (paired-end)."
    echo "  -r  SINGLE_FQ      Path to the single-end FASTQ file."
    exit 1
}

# --- Main Function ---
main() {
    # --- Parse Command-Line Arguments ---
    READ1=""
    READ2=""
    SINGLE_READ=""
    SAMPLE_ID=""

    while getopts "s:1:2:r:" opt; do
        case ${opt} in
            s ) SAMPLE_ID=$OPTARG ;;
            1 ) READ1=$OPTARG ;;
            2 ) READ2=$OPTARG ;;
            r ) SINGLE_READ=$OPTARG ;;
            \? ) usage ;;
        esac
    done

    # Validate arguments
    if [ -z "$SAMPLE_ID" ]; then
        log "Error: Sample ID is required." >&2; usage
    fi
    if [ -z "$READ1" ] && [ -z "$SINGLE_READ" ]; then
        log "Error: You must provide either paired-end reads (-1 and -2) or a single-end read (-r)." >&2; usage
    fi
    if [ -n "$READ1" ] && [ -z "$READ2" ]; then
        log "Error: Paired-end requires both -1 and -2 files." >&2; usage
    fi

    # --- Configuration Section ---
    # Customize these paths for your system.
    
    # Databases
    HOST_GENOME="/path/to/db/hg38/hg38.fa.gz"
    DVF_MODELS="/path/to/db/DeepVirFinder/models"
    GENOMAD_DB="/path/to/db/genomad_db"
    CHECKV_DB="/path/to/db/checkv-db-v1.4"
    
    # Parameters
    THREADS=8
    MEMORY_GB=50 # Memory for SPAdes in GB
    MIN_CONTIG_LEN=10000 # Minimum contig length for viral identification

    # Main Output Directory (created relative to script execution location)
    OUT_DIR="./${SAMPLE_ID}_analysis"

    # --- Tool Definitions ---
    # Ensure these tools are in your PATH or provide the full path.
    local FASTP="fastp"
    local SOAPNUKE="SOAPnuke"
    local BWA="bwa"
    local SAMTOOLS="samtools"
    local BAMTOFASTQ="bamToFastq"
    local MEGAHIT="megahit"
    local SEQKIT="seqkit"
    local DVF="dvf"
    local DVF_EXTRACT="dvf.extract.py"
    local VIRSORTER="virsorter"
    local GENOMAD="genomad"
    local CDHIT="cd-hit-est"
    local CHECKV="checkv"
    
    # --- Start Pipeline ---
    log "--- Viro-Flow Stage 1 Initialized for Sample: $SAMPLE_ID ---"
    mkdir -p "$OUT_DIR"

    # --- 1. Quality Control ---
    run_qc

    # --- 2. Host Read Removal ---
    remove_host

    # --- 3. Assembly & Filtering ---
    assemble_and_filter

    # --- 4. Viral Identification & QC ---
    identify_viruses

    log "--- Viro-Flow Stage 1 Finished for Sample: $SAMPLE_ID ---"
    log "Final high-quality viral sequences are in: ${OUT_DIR}/${SAMPLE_ID}_final_viruses.fna"
}

# --- Pipeline Functions ---

run_qc() {
    log "[1/4] Running Quality Control..."
    local QC_DIR="${OUT_DIR}/01_qc"
    local FASTP_DIR="${QC_DIR}/fastp"
    local SOAPNUKE_DIR="${QC_DIR}/soapnuke"
    mkdir -p "$FASTP_DIR" "$SOAPNUKE_DIR"

    if [ -n "$READ1" ]; then # Paired-end
        $FASTP -i "$READ1" -I "$READ2" \
               -o "${FASTP_DIR}/fastp_1.fq.gz" -O "${FASTP_DIR}/fastp_2.fq.gz" \
               -5 -3 -z 4 -q 20 -c -l 30 \
               -j "${QC_DIR}/${SAMPLE_ID}.fastp.json" -h "${QC_DIR}/${SAMPLE_ID}.fastp.html"

        $SOAPNUKE filter -1 "${FASTP_DIR}/fastp_1.fq.gz" -2 "${FASTP_DIR}/fastp_2.fq.gz" \
                         -C "${SOAPNUKE_DIR}/soapnuke_1.fq.gz" -D "${SOAPNUKE_DIR}/soapnuke_2.fq.gz" -o "$SOAPNUKE_DIR"
        
        # Set global variables for next steps
        CLEAN_READ1="${SOAPNUKE_DIR}/soapnuke_1.fq.gz"
        CLEAN_READ2="${SOAPNUKE_DIR}/soapnuke_2.fq.gz"
    else # Single-end
        $FASTP -i "$SINGLE_READ" -o "${FASTP_DIR}/fastp.fq.gz" \
               -5 -3 -z 4 -q 20 -c -l 30 \
               -j "${QC_DIR}/${SAMPLE_ID}.fastp.json" -h "${QC_DIR}/${SAMPLE_ID}.fastp.html"

        $SOAPNUKE filter -f "${FASTP_DIR}/fastp.fq.gz" -C "${SOAPNUKE_DIR}/soapnuke.fq.gz" -o "$SOAPNUKE_DIR"
        
        CLEAN_READ1="${SOAPNUKE_DIR}/soapnuke.fq.gz"
        CLEAN_READ2="" # Ensure READ2 is empty for single-end
    fi
    log "QC complete. Clean reads are in ${SOAPNUKE_DIR}"
}

remove_host() {
    log "[2/4] Removing Host Genome..."
    local RMHOST_DIR="${OUT_DIR}/02_rmhost"
    mkdir -p "$RMHOST_DIR"
    
    if [ -n "$CLEAN_READ2" ]; then # Paired-end
        $BWA mem -t "$THREADS" "$HOST_GENOME" "$CLEAN_READ1" "$CLEAN_READ2" | \
        $SAMTOOLS view -b -f 12 -F 256 | \
        $SAMTOOLS sort -n -@ "$THREADS" | \
        $BAMTOFASTQ -i /dev/stdin -fq "${RMHOST_DIR}/rmhost_1.fastq.gz" -fq2 "${RMHOST_DIR}/rmhost_2.fastq.gz"

        # Set global variables
        RMHOST_READ1="${RMHOST_DIR}/rmhost_1.fastq.gz"
        RMHOST_READ2="${RMHOST_DIR}/rmhost_2.fastq.gz"
    else # Single-end
        $BWA mem -t "$THREADS" "$HOST_GENOME" "$CLEAN_READ1" | \
        $SAMTOOLS view -b -f 4 -F 256 | \
        $BAMTOFASTQ -i /dev/stdin -fq "${RMHOST_DIR}/rmhost.fastq.gz"

        RMHOST_READ1="${RMHOST_DIR}/rmhost.fastq.gz"
        RMHOST_READ2=""
    fi
    log "Host removal complete. Non-host reads are in ${RMHOST_DIR}"
}

assemble_and_filter() {
    log "[3/4] Assembling and Filtering Contigs..."
    local ASSEMBLY_DIR="${OUT_DIR}/03_assembly"
    mkdir -p "$ASSEMBLY_DIR"

    if [ -n "$RMHOST_READ2" ]; then # Paired-end
        $MEGAHIT -1 "$RMHOST_READ1" -2 "$RMHOST_READ2" -t "$THREADS" --presets meta-large -o "$ASSEMBLY_DIR"
    else # Single-end
        $MEGAHIT -r "$RMHOST_READ1" -t "$THREADS" --presets meta-large -o "$ASSEMBLY_DIR"
    fi
    mv "${ASSEMBLY_DIR}/final.contigs.fa" "${ASSEMBLY_DIR}/assembly.fa"
    
    # Filter
    $SEQKIT seq -m "$MIN_CONTIG_LEN" "${ASSEMBLY_DIR}/assembly.fa" > "${ASSEMBLY_DIR}/contigs_filtered.fa"
    $SEQKIT stats -j "$THREADS" -a "${ASSEMBLY_DIR}/assembly.fa" > "${ASSEMBLY_DIR}/assembly_stats.txt"
    log "Assembly complete. Filtered contigs are in ${ASSEMBLY_DIR}"
}

identify_viruses() {
    log "[4/4] Identifying and Filtering Viral Sequences..."
    local VIRAL_DIR="${OUT_DIR}/04_viral_id"
    local DVF_OUT="${VIRAL_DIR}/dvf"
    local VS2_OUT="${VIRAL_DIR}/vs2"
    local GENOMAD_OUT="${VIRAL_DIR}/genomad"
    local CHECKV_DIR="${OUT_DIR}/05_checkv"
    mkdir -p "$DVF_OUT" "$VS2_OUT" "$GENOMAD_OUT" "$CHECKV_DIR"
    local FILTERED_CONTIGS="${OUT_DIR}/03_assembly/contigs_filtered.fa"

    log "  -> Running DeepVirFinder..."
    $DVF -i "$FILTERED_CONTIGS" -o "$DVF_OUT" -m "$DVF_MODELS" -l "$MIN_CONTIG_LEN" -c "$THREADS"
    $DVF_EXTRACT "${DVF_OUT}/*gt*bp_dvfpred.txt" "${DVF_OUT}/dvfpred.id"
    $SEQKIT grep -n -f "${DVF_OUT}/dvfpred.id" "$FILTERED_CONTIGS" > "${DVF_OUT}/dvfpred.fa"

    log "  -> Running VirSorter2..."
    $VIRSORTER run -w "$VS2_OUT" -i "$FILTERED_CONTIGS" --min-length "$MIN_CONTIG_LEN" -j "$THREADS" all

    log "  -> Running GeNoMad..."
    $GENOMAD end-to-end --cleanup --splits "$THREADS" "$FILTERED_CONTIGS" "$GENOMAD_OUT" "$GENOMAD_DB"

    local COMBINED_FA="${VIRAL_DIR}/combined_viral.fa"
    local CDHIT_FA="${VIRAL_DIR}/combined_viral_cdhit99.fa"
    cat "${VS2_OUT}/final-viral-combined.fa" "${DVF_OUT}/dvfpred.fa" "${GENOMAD_OUT}"/*_summary/*virus.fna > "$COMBINED_FA"
    $CDHIT -i "$COMBINED_FA" -o "$CDHIT_FA" -c 0.99 -aL 0.9 -M 0 -T "$THREADS"

    log "  -> Running CheckV for quality assessment..."
    $CHECKV end_to_end "$CDHIT_FA" "$CHECKV_DIR" -d "$CHECKV_DB" -t "$THREADS"

    awk '($10 >= 50) || ($8 == "Not-determined" && $2 >= 30000)' "${CHECKV_DIR}/quality_summary.tsv" | \
    grep -v "contig_id" | cut -f1 > "${CHECKV_DIR}/quality.id"
    $SEQKIT grep -n -f "${CHECKV_DIR}/quality.id" "$CDHIT_FA" > "${CHECKV_DIR}/high_quality_viral.fa"

    local FINAL_VIRAL_FASTA="${OUT_DIR}/${SAMPLE_ID}_final_viruses.fna"
    sed "s/^>/>${SAMPLE_ID}_/" "${CHECKV_DIR}/high_quality_viral.fa" > "$FINAL_VIRAL_FASTA"
    log "Viral identification complete."
}


# --- Execute Main Function ---
main "$@"
