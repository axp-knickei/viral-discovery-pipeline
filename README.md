# **Viro-Flow (viral-discovery-pipeline): A Viral Metagenomics Discovery Pipeline**

Viro-Flow is a robust and automated bioinformatics pipeline for the identification and analysis of viral sequences from raw metagenomic sequencing data. It processes one or more samples from raw reads to a final abundance matrix of viral operational taxonomic units (vOTUs), making it suitable for comparative viral metagenomics.

[Image of a DNA sequencing workflow diagram](https://encrypted-tbn0.gstatic.com/licensed-image?q=tbn:ANd9GcQKvDI7XuUGWnKOtnd5K7glwtGeCbWnQHj-ufKkR2i6qXj_XW8Xu0gIvCnOdskoRk0w4jJRTf3YHvPHAMzy223wmUzV5zhQDh-xPUFIVkF2DEV7I6k)

The pipeline is structured into two main stages:

1. **viral\_pipeline.sh**: A per-sample processing script that takes raw FASTQ files and performs QC, host read removal, assembly, and viral contig identification.  
2. **combine\_and\_analyze.sh**: A multi-sample script that aggregates the results from all samples, performs species-level clustering to define vOTUs, and calculates their relative abundance across all samples.

## **Features**

* **Automated & Robust**: The pipeline is designed to run from start to finish with minimal intervention. It includes error checking to halt on command failures.  
* **Configurable**: All tool paths, databases, and key parameters are located in a central configuration section for easy setup.  
* **Scalable**: Efficiently processes individual samples and then aggregates them for large-scale comparative analysis.  
* **Best Practices**: Utilizes a suite of well-regarded bioinformatics tools for each step of the analysis.  
* **Detailed Documentation**: Includes comprehensive information on dependencies and tool functions.

## **Installation & Dependencies**

To use Viro-Flow, you must have several bioinformatics tools installed and available in your system's PATH. We recommend using a Conda environment to manage these dependencies.

### **1\. Create a Conda Environment**

conda create \-n viroflow-env \-c conda-forge \-c bioconda \\  
    fastp soapnuke bwa samtools bedtools \\  
    spades megahit seqkit python=3.8 \\  
    virsorter2 'dvf-py=1.0' genomad checkv \\  
    fastani bowtie2 coverm

*Note: Some tools like bamToFastq may require separate installation if not available via Conda.*

### **2\. Activate the Environment**

conda activate viroflow-env

### **3\. Download Databases**

You will also need to download the required databases and specify their paths in the configuration section of the scripts. See the [Tools Documentation](https://www.google.com/search?q=./DOCS/TOOLS.md) for links and details.

## **Usage**

### **Stage 1: Per-Sample Processing**

Run viral\_pipeline.sh for each of your samples. This script will create a dedicated output folder for each sample.

**For paired-end data:**

./viral\_pipeline.sh \\  
    \-s MySample1 \\  
    \-1 /path/to/reads\_1.fastq.gz \\  
    \-2 /path/to/reads\_2.fastq.gz

**For single-end data:**

./viral\_pipeline.sh \\  
    \-s MySample2 \\  
    \-r /path/to/reads.fastq.gz

After this stage, each sample directory will contain a \*\_final\_viruses.fna file with high-quality viral contigs.

### **Stage 2: Multi-Sample Analysis**

Once all samples have been processed, run combine\_and\_analyze.sh to cluster the viruses and calculate their abundance.

./combine\_and\_analyze.sh \\  
    \-i /path/to/parent/dir/of/all/sample/folders \\  
    \-o /path/to/final\_analysis\_output

### **Key Outputs**

* **vOTUs\_representatives.fa**: A FASTA file containing the representative sequence for each viral "species" (vOTU). This is your final viral database.  
* **vOTU\_abundance\_matrix\_tpm.tsv**: A tab-separated file with vOTUs as rows and samples as columns, containing the TPM (Transcripts Per Million) abundance values. This table is ready for downstream statistical analysis and visualization.

## **Citation**

If you use this pipeline, please cite the individual tools used in the workflow. You can find a complete list and links in the [Tools Documentation](https://www.google.com/search?q=./DOCS/TOOLS.md).