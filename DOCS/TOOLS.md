# **Viro-Flow: Tool and Dependency Documentation**

This document provides a detailed overview of the external bioinformatics tools used in the Viro-Flow pipeline, along with their functions and links to their official sources.

## **Stage 1: viral\_pipeline.sh**

### **Data Quality Control**

* **fastp**  
  * **Function**: An ultra-fast all-in-one FASTQ preprocessor for quality control, adapter trimming, and filtering.  
  * **Key Parameters**: \-q 20 (minimum quality score), \-l 30 (minimum read length).  
  * **Source**: [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp)  
* **SOAPnuke**  
  * **Function**: A tool for filtering sequencing reads, used here as a secondary QC step.  
  * **Source**: [https://github.com/BGI-flexlab/SOAPnuke](https://github.com/BGI-flexlab/SOAPnuke)

### **Host Read Removal**

* **BWA (Burrows-Wheeler Aligner)**  
  * **Function**: A fast and accurate short-read aligner. Used to map reads to a host genome (e.g., human hg38) for identification and removal.  
  * **Key Parameters**: mem (algorithm for longer reads), \-t (threads).  
  * **Source**: [https://github.com/lh3/bwa](https://github.com/lh3/bwa)  
* **SAMtools**  
  * **Function**: A suite of utilities for interacting with high-throughput sequencing data alignments in SAM/BAM format. Used here to filter for unmapped reads.  
  * **Key Parameters**: view \-f 12 (select properly paired, unmapped reads), sort (sort alignments).  
  * **Source**: [https://github.com/samtools/samtools](https://github.com/samtools/samtools)  
* **bamToFastq**  
  * **Function**: Converts a BAM file back to FASTQ format.  
  * **Source**: Part of the bedtools suite: [https://github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2)

### **Assembly**

* **MEGAHIT**  
  * **Function**: A memory-efficient and fast metagenomic assembler. It reconstructs long contigs from short reads.  
  * **Key Parameters**: \--presets meta-large (optimizes for complex metagenomes).  
  * **Source**: [https://github.com/voutcn/megahit](https://github.com/voutcn/megahit)

### **Viral Identification & QC**

* **SeqKit**  
  * **Function**: A versatile toolkit for FASTA/Q file manipulation. Used for filtering contigs by length and extracting sequences.  
  * **Key Parameters**: seq \-m (filter by minimum length), grep (extract sequences by ID).  
  * **Source**: [https://github.com/shenwei356/seqkit](https://github.com/shenwei356/seqkit)  
* **DeepVirFinder (DVF)**  
  * **Function**: A deep learning-based tool for predicting whether a DNA sequence is of viral origin.  
  * **Source**: [https://github.com/jessieren/DeepVirFinder](https://github.com/jessieren/DeepVirFinder)  
* **VirSorter2**  
  * **Function**: A widely-used tool that identifies viral sequences in metagenomic data using a guilt-by-association approach with viral protein families.  
  * **Source**: [https://github.com/jiarong/VirSorter2](https://github.com/jiarong/VirSorter2)  
* **GeNoMad**  
  * **Function**: Identifies viruses and other mobile genetic elements in genomic or metagenomic data.  
  * **Source**: [https://github.com/apcamargo/genomad](https://github.com/apcamargo/genomad)  
* **CD-HIT**  
  * **Function**: A program for clustering and comparing protein or nucleotide sequences. Used here to dereplicate sequences at 99% identity.  
  * **Key Parameters**: \-c 0.99 (sequence identity threshold).  
  * **Source**: [https://github.com/weizhongli/cdhit](https://github.com/weizhongli/cdhit)  
* **CheckV**  
  * **Function**: Assesses the quality of viral contigs by estimating genome completeness and identifying potential host contamination.  
  * **Source**: [https://bitbucket.org/berkeleylab/checkv](https://bitbucket.org/berkeleylab/checkv)

## **Stage 2: combine\_and\_analyze.sh**

### **Clustering**

* **fastANI**  
  * **Function**: A fast alignment-free tool for computing whole-genome Average Nucleotide Identity (ANI). It is much faster than BLAST-based methods for species-level clustering.  
  * **Source**: [https://github.com/ParBLiSS/FastANI](https://github.com/ParBLiSS/FastANI)

### **Abundance Calculation**

* **Bowtie2**  
  * **Function**: A fast and memory-efficient tool for aligning sequencing reads to long reference sequences. Used here to map reads to the final vOTU database.  
  * **Key Parameters**: \--very-sensitive (preset for high accuracy).  
  * **Source**: [https://github.com/BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2)  
* **CoverM**  
  * **Function**: A tool for calculating coverage of genomes or contigs from read alignments. It can calculate various normalized abundance metrics.  
  * **Key Parameters**: \-m tpm (calculates Transcripts Per Million).  
  * **Source**: [https://github.com/wwood/CoverM](https://github.com/wwood/CoverM)