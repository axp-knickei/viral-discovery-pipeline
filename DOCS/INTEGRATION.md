# Viro-Flow & Cataloger: The Complete Viral Discovery Workflow

This document details the relationship between the **Viral Discovery Pipeline (Viro-Flow)** and the **Viral Genome Cataloger**, and how to integrate them for a complete viral discovery and cataloging workflow.

## Executive Summary

* **viral-discovery-pipeline (Viro-Flow)** is a comprehensive, end-to-end **metagenomics pipeline**. It takes raw sequencing data (FASTQ) and performs every step necessary to find viruses, including cleaning, assembly, identification, and abundance profiling.
* **viral-genome-cataloger** is a specialized **clustering tool**. It takes already-assembled viral genomes (FASTA) and groups them into non-redundant species catalogs (vOTUs) using a specific algorithm to prevent clustering artifacts.

## 1. Purpose and Scope

| Feature | Viral Discovery Pipeline (Viro-Flow) | Viral Genome Cataloger |
| --- | --- | --- |
| **Primary Goal** | **Discovery**: Turn raw sequencing reads into a viral community abundance matrix. | **Dereplication**: Turn a collection of viral genomes into a non-redundant catalog. |
| **Input** | Raw sequencing data (**FASTQ** files). | Assembled genome sequences (**FASTA** files). |
| **Output** | Quality reports, viral contigs, and a sample-by-virus abundance matrix (TPM). | A single FASTA file of representative viral species (vOTUs) and a cluster map. |
| **Use Case** | You have raw data from a sequencer and want to know "What viruses are in my samples?" | You have viral genomes from multiple studies/samples and want to merge them into a single database. |

## 2. Methodological Differences & Scientific Rationale

The most significant difference lies in how they handle **Clustering** (grouping similar viruses into species).

### Viro-Flow: Connected Components (The "Miner")
* **Method:** Uses a "connected components" approach.
* **Tools:** Uses **CD-HIT-EST** for initial dereplication (99% identity) followed by **FastANI** for species-level clustering.
* **Algorithm:** It calculates Average Nucleotide Identity (ANI). If Genome A matches Genome B, they are linked. It uses `networkx` to find connected components in this graph.
* **Pros/Cons:** Standard approach, but susceptible to "chaining" (where A matches B, and B matches C, so A and C are grouped even if they are not similar).

### Genome Cataloger: Greedy Star-Topology (The "Librarian")
* **Method:** Uses a "Greedy Star-Topology" approach.
* **Tools:** Uses **skani** (a newer, faster alternative to FastANI).
* **Algorithm:**
    1. Sorts all genomes by length.
    2. Picks the longest as a "centroid."
    3. Recruits only genomes that match the centroid directly.
    4. Removes them from the pool and repeats.
* **Scientific Rationale:** This method is explicitly designed to prevent "chaining" artifacts, ensuring all members of a cluster are similar to the representative. It provides a stricter, scientifically rigorous species catalog.

## 3. Workflow and Architecture

**Viral Discovery Pipeline (Viro-Flow):**
* **Complexity:** High. Orchestrates a massive suite of bioinformatics tools including fastp, BWA, MEGAHIT, VirSorter2, GeNoMad, and CheckV.
* **Implementation:** Hybrid (Bash scripts and Python CLI).
* **Environment:** Reliance on Conda/Docker due to complex dependencies.

**Viral Genome Cataloger:**
* **Complexity:** Low. Focuses on one task efficiently.
* **Implementation:** Pure Python package wrapping external binaries.
* **Dependencies:** Lightweight (skani, SeqKit).

## 4. Integration Guide: From Discovery to Catalog

The combined workflow follows a "Discovery then Dereplication" logic.

### Step 1: Mining with Viro-Flow
Run the `viral_pipeline.sh` for each of your sequencing samples. This stage handles quality control, host removal, assembly, and viral identification.

* **Primary Output**: `SampleID_final_viruses.fna`.
* **Location**: Found in the individual analysis directory for each sample.

### Step 2: The Handover
Collect the `.fna` files from all your Viro-Flow analysis directories and place them into a single input folder for the cataloger.

```bash
# Example of gathering outputs
mkdir all_discovered_viruses
cp *_analysis/*_final_viruses.fna all_discovered_viruses/
```

### Step 3: Cataloging with Genome Cataloger
Run the `viral-cataloger` on your gathered sequences. It uses `skani` for high-resolution ANI estimation and a Greedy Star-Topology algorithm.

* **Scientific Standard**: Defaults to **95% ANI** and **85% Coverage** (international standard).
* **Primary Output**: `catalog_vOTU_catalog.fasta` (one representative per species).

## Why use the Cataloger instead of Viro-Flow's built-in clustering?

Viro-Flow includes `combine_and_analyze.sh` which performs clustering. However, `viral-genome-cataloger` is recommended for the final reference database because:

1.  **Speed:** `skani` is generally faster than the **FastANI** method used in the main pipeline.
2.  **Accuracy:** "Star-topology" clustering is stricter and prevents "chaining" errors (where dissimilar viruses get grouped together through a chain of intermediates).
3.  **Flexibility:** You can easily combine your results with viral genomes found in *other* studies or public databases. You can just dump those FASTA files into the input folder along with your pipeline results, and `viral-genome-cataloger` will process them all together.

---

**Recommendation:**
* Use **`viral-discovery-pipeline`** for the heavy lifting of assembly and viral identification from raw data.
* Use **`viral-genome-cataloger`** if you have already run the discovery pipeline on multiple different datasets and now need to combine all your results into a single, high-quality reference database without duplicates.
