# HPC User Guide for Viro-Flow

This guide explains how to run the **Viro-Flow** pipeline on High Performance Computing (HPC) clusters. HPC environments typically use job schedulers (like Slurm or PBS) and require specific methods for software installation (Containers or Conda).

## Option 1: Using Singularity / Apptainer (Recommended)

Containers are the best way to ensure reproducibility and avoid installation headaches on shared clusters.

### 1. Build the Image
You usually cannot build images directly on the HPC (requires root). Build it on your local machine or use a remote builder, then transfer it.

**On your local machine (with Singularity/Apptainer installed):**
```bash
sudo singularity build viroflow.sif Singularity.def
```

**Transfer to HPC:**
```bash
scp viroflow.sif user@hpc.institute.edu:/path/to/project/
```

### 2. Run an Interactive Shell (Testing)
```bash
singularity shell viroflow.sif
# You are now inside the container
viral-pipeline --help
exit
```

### 3. Submit a Slurm Job
Create a file named `submit_viroflow.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=viroflow
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Load Singularity module if required
module load singularity

# Set paths
SIF_IMAGE="/path/to/project/viroflow.sif"
INPUT_DIR="/path/to/project/data"
OUTPUT_DIR="/path/to/project/results"
DB_DIR="/path/to/project/databases"

# Create output directories
mkdir -p logs

# Run the pipeline
# We bind mount the data directories so the container can see them
singularity exec \
    --bind $INPUT_DIR:/data \
    --bind $OUTPUT_DIR:/outputs \
    --bind $DB_DIR:/databases \
    $SIF_IMAGE \
    viral-pipeline run \
    --input-dir /data \
    --output-dir /outputs \
    --config-path config.yaml
```

Submit it:
```bash
sbatch submit_viroflow.sh
```

---

## Option 2: Using Conda / Mamba

If you cannot use containers, you can install everything into a Conda environment in your home or project directory.

### 1. Create Environment
```bash
# Load conda/mamba module
module load miniconda3

# Create environment with all dependencies
mamba create -n viroflow -c conda-forge -c bioconda \
    python=3.11 \
    bowtie2 coverm cd-hit fastani samtools bedtools \
    fastp bwa megahit virsorter=2 checkv genomad seqkit \
    pandas networkx

# Activate
source activate viroflow

# Install the pipeline package
pip install .
```

### 2. Submit a Slurm Job
Create `submit_viroflow_conda.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=viroflow_conda
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Load modules and activate environment
module load miniconda3
source activate viroflow

# Run
viral-pipeline run \
    --input-dir data/ \
    --output-dir results/ \
    --config-path config.yaml
```

---

## Performance Tips for HPC

1.  **Threads**: The pipeline is parallelized. Match `--threads` (in config) to `#SBATCH --cpus-per-task`.
    *   *Config*: `tools: threads: 16`
    *   *Slurm*: `#SBATCH --cpus-per-task=16`
2.  **Memory**: Viral assembly (MEGAHIT) and clustering can be memory intensive.
    *   Start with **32GB - 64GB** for standard metagenomes.
    *   Use **100GB+** for deep sequencing or soil samples.
3.  **IO**: Read/Write operations can be slow on network storage.
    *   If available, use fast scratch storage (e.g., `/scratch` or `$TMPDIR`) for the `output_dir`, then copy results to permanent storage at the end of the script.

```bash
# Example using scratch
export OUT_DIR=$TMPDIR/viroflow_results

viral-pipeline run ... --output-dir $OUT_DIR

# Copy back
cp -r $OUT_DIR /path/to/project/final_results/
```
