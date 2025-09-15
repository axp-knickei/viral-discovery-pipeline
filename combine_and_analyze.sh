#!/bin/bash

# Viro-Flow Stage 2: Multi-Sample Clustering and Abundance Analysis
#
# This script aggregates viral contigs from multiple samples, clusters them
# into vOTUs (species), and calculates their abundance across all samples.
# Run this script after processing all samples with 'viral_pipeline.sh'.

# --- Script Setup ---
set -e
set -u
set -o pipefail

# --- Logging Function ---
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] -" "$1"
}

# --- Usage Function ---
usage() {
    echo "Usage: $0 -i <INPUT_DIR> -o <OUTPUT_DIR>"
    echo "  -i  INPUT_DIR    Path to the directory containing all individual sample analysis folders."
    echo "  -o  OUTPUT_DIR   Path to the directory where combined analysis results will be saved."
    exit 1
}

# --- Main Function ---
main() {
    # --- Parse Arguments ---
    INPUT_DIR=""
    OUTPUT_DIR=""

    while getopts "i:o:" opt; do
        case ${opt} in
            i ) INPUT_DIR=$OPTARG ;;
            o ) OUTPUT_DIR=$OPTARG ;;
            \? ) usage ;;
        esac
    done

    if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
        log "Error: Both input and output directories are required." >&2; usage
    fi

    # --- Configuration ---
    local THREADS=16
    local ANI_THRESHOLD=0.95 # 95% ANI for species clustering
    local MIN_COVERAGE=0.85  # 85% coverage for clustering

    # --- Tool Definitions ---
    local SEQKIT="seqkit"
    local CDHIT="cd-hit-est"
    local FASTANI="fastANI"
    local BOWTIE2_BUILD="bowtie2-build"
    local BOWTIE2="bowtie2"
    local SAMTOOLS="samtools"
    local COVERM="coverm"

    # --- Start Pipeline ---
    log "--- Viro-Flow Stage 2 Initialized ---"
    mkdir -p "$OUTPUT_DIR"

    # --- 1. Aggregate and Dereplicate ---
    aggregate_and_dereplicate

    # --- 2. Cluster into vOTUs ---
    cluster_vot_us

    # --- 3. Calculate Abundance ---
    calculate_abundance

    log "--- Viro-Flow Stage 2 Finished! ---"
    log "Key output files are located in: ${OUTPUT_DIR}"
}

# --- Pipeline Functions ---

aggregate_and_dereplicate() {
    log "[1/3] Aggregating and dereplicating all viral sequences..."
    local AGGREGATE_DIR="${OUTPUT_DIR}/01_aggregation"
    mkdir -p "$AGGREGATE_DIR"
    ALL_VIRUSES_FA="${AGGREGATE_DIR}/all_samples_viruses.fna"
    DEREPLICATED_FA="${AGGREGATE_DIR}/all_viruses_dereplicated.fna"

    find "${INPUT_DIR}" -name "*_final_viruses.fna" -print0 | xargs -0 cat > "$ALL_VIRUSES_FA"
    $CDHIT -i "$ALL_VIRUSES_FA" -o "$DEREPLICATED_FA" -c 0.99 -n 10 -M 0 -T "$THREADS"
    log "Aggregation complete. Dereplicated sequences at ${DEREPLICATED_FA}"
}

cluster_vot_us() {
    log "[2/3] Performing species clustering (vOTUs)..."
    local CLUSTERING_DIR="${OUTPUT_DIR}/02_clustering"
    mkdir -p "$CLUSTERING_DIR"
    local FASTANI_OUT="${CLUSTERING_DIR}/fastani_results.txt"
    local CLUSTER_MAP="${CLUSTERING_DIR}/votu_cluster_map.tsv"
    VOTU_REPS_FA="${CLUSTERING_DIR}/vOTUs_representatives.fa"

    $FASTANI --ql "$DEREPLICATED_FA" --rl "$DEREPLICATED_FA" -o "$FASTANI_OUT" -t "$THREADS"

    awk -v ani="$ANI_THRESHOLD" -v cov="$MIN_COVERAGE" '$3/100 >= ani && $4/$5 >= cov {print $1 "\t" $2}' "$FASTANI_OUT" | \
    python3 -c 'import networkx as nx, sys; G=nx.from_pandas_edgelist(pd.read_csv(sys.stdin, sep="\t", header=None), 0, 1); [print(f"{node}\tvOTU_{i+1}") for i, cluster in enumerate(nx.connected_components(G)) for node in cluster]' > "$CLUSTER_MAP"

    $SEQKIT grep -f <(cut -f1 "$CLUSTER_MAP") "$DEREPLICATED_FA" | \
    $SEQKIT fx2tab -nl | sort -k1,1 | join -1 1 -2 1 <(sort -k1,1 "$CLUSTER_MAP") -o 2.2,1.1,1.2 | \
    sort -k1,1 -k3,3nr | sort -u -k1,1 --merge | cut -f2 | \
    $SEQKIT grep -f - "$DEREPLICATED_FA" > "$VOTU_REPS_FA"
    
    log "Clustering complete. Found $(grep -c '>' "$VOTU_REPS_FA") vOTUs."
}

calculate_abundance() {
    log "[3/3] Calculating relative abundance..."
    local INDEX_DIR="${OUTPUT_DIR}/03_abundance_index"
    local ABUNDANCE_DIR="${OUTPUT_DIR}/04_abundance"
    local TMP_DIR="${ABUNDANCE_DIR}/per_sample_tpm"
    mkdir -p "$INDEX_DIR" "$TMP_DIR"

    $BOWTIE2_BUILD "$VOTU_REPS_FA" "${INDEX_DIR}/vOTUs_db"

    for SAMPLE_DIR in "${INPUT_DIR}"/*_analysis; do
        if [ -d "$SAMPLE_DIR" ]; then
            local SAMPLE_ID=$(basename "$SAMPLE_DIR" | sed 's/_analysis$//')
            log "  -> Mapping reads for sample: $SAMPLE_ID"
            local RMHOST_DIR="${SAMPLE_DIR}/02_rmhost"
            local READ1="${RMHOST_DIR}/rmhost_1.fastq.gz"
            local READ2="${RMHOST_DIR}/rmhost_2.fastq.gz"
            local BAM_OUT="${TMP_DIR}/${SAMPLE_ID}.bam"

            if [ -f "$READ2" ]; then # Paired-end
                $BOWTIE2 -p "$THREADS" -x "${INDEX_DIR}/vOTUs_db" -1 "$READ1" -2 "$READ2" --very-sensitive | \
                $SAMTOOLS view -bS | $SAMTOOLS sort -@ "$THREADS" -o "$BAM_OUT"
            else # Single-end
                 $BOWTIE2 -p "$THREADS" -x "${INDEX_DIR}/vOTUs_db" -U "$READ1" --very-sensitive | \
                $SAMTOOLS view -bS | $SAMTOOLS sort -@ "$THREADS" -o "$BAM_OUT"
            fi
            $COVERM contig -b "$BAM_OUT" -m tpm --min-covered-fraction 0.60 > "${TMP_DIR}/${SAMPLE_ID}.tpm.txt"
        fi
    done

    log "  -> Aggregating abundance results into a matrix..."
    python3 -c 'import pandas as pd, sys; all_files = sys.argv[1:]; dfs = [pd.read_csv(f, sep="\t", header=0, index_col=0, names=[f.split("/")[-1].replace(".tpm.txt","")]) for f in all_files]; df = pd.concat(dfs, axis=1).fillna(0); df.index.name = "vOTU"; df.to_csv(sys.stdout, sep="\t")' "${TMP_DIR}"/*.tpm.txt > "${ABUNDANCE_DIR}/vOTU_abundance_matrix_tpm.tsv"
    
    log "Abundance calculation complete."
}

# --- Execute Main Function ---
main "$@"
