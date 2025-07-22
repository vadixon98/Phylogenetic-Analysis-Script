# AOS Phylogeny & Functional Analysis Pipeline

**Script:** `Final_Complete_Annotated.R`
**Author:** Victoria A. Dixon
**Date:** 2025-05-02

## Description

This R pipeline automates retrieval and analysis of Allene Oxide Synthase (AOS; CYP74A) genes in *Solanum lycopersicum* (tomato). It performs:

1. **Sequence Retrieval**: Downloads coding and protein sequences, filters for AOS.
2. **Promoter Extraction**: Fetches \~2 kb upstream sequences for AOS genes.
3. **Multiple Sequence Alignment & Motif Scanning**: Aligns protein sequences and scans for motifs (e.g., heme-binding, jasmonate elements).
4. **Phylogenetic Reconstruction**: Builds neighbor‑joining and maximum‑likelihood trees with bootstrap support, and saves publication‑quality figures.

All intermediate and final outputs are organized under `data/` and `results/` directories.

---

## Requirements

* **R** ≥ 4.0
* Internet connection for sequence and annotation downloads
* External tools installed if using MEME/FIMO (optional for motif scanning)

No non‑standard R packages are required; dependencies are installed automatically via CRAN and Bioconductor.

---

## Installation & Setup

1. Clone or download the repository:

   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```
2. Open R or RStudio in this directory.
3. Run the pipeline script:

   ```r
   source("Final_Complete_Annotated.R")
   ```

   This will:

   * Create `data/` and `results/` if missing
   * Install and load required CRAN/Bioconductor packages
   * Execute all pipeline steps in sequence

---

## Pipeline Steps

### 1. Directory Setup

* Creates `data/` and `results/` for organized output storage.

### 2. Package Management

* Installs missing CRAN packages (`biomartr`, `ape`, `phangorn`, `tidyverse`, `seqinr`).
* Installs missing Bioconductor packages (`rentrez`, `msa`, `Biostrings`, `ggtree`, `biomaRt`, `rtracklayer`, `IRanges`).
* Loads all dependencies.

### 3. Sequence Retrieval

* Downloads CDS and proteome for tomato using **biomartr** and **rentrez**.
* Filters for AOS/CYP74A sequences by header patterns.
* Saves FASTA files: `data/Sl_AOS_CDS.fasta` and `data/Sl_AOS_protein.fasta`.

### 4. Promoter Extraction

* Imports GFF3 annotation and genome FASTA.
* Identifies AOS gene features.
* Extracts \~2 kb upstream promoter regions (reverse-complemented for negative strand).
* Saves to `data/AOS_promoters_2kb.fasta`.

### 5. Alignment & Motif Scanning

* Aligns AOS protein sequences with **msa::msa** (ClustalW).
* Exports alignment: `results/AOS_protein_aln.fasta`.
* Scans for heme-binding motif (`F.GGPRC`) and prints hit positions.
* Runs external FIMO motif scan on promoters (requires MEME suite) and outputs to `results/fimo_JA_elements/`.

### 6. Phylogenetic Reconstruction

* Converts MSA to **phangorn** phyDat object.
* Computes distance matrix and neighbor-joining tree.
* Performs maximum-likelihood optimization (JTT+Γ+I).
* Bootstraps ML tree (1000 replicates) and annotates node support.
* Saves trees as PNG:

  * `results/AOS_NJ_tree.png`
  * `results/AOS_ML_tree.png`

### 7. Pipeline Completion

* Final message indicating successful run and location of outputs.

---

## Usage Tips

* **Re-running**: To update with new data, delete contents of `data/` and `results/` before rerunning.
* **Parameter Tweaks**: Modify gap penalties, bootstrap replicates, or motif files by editing the script variables.
* **Headless Runs**: Use `Rscript Final_Complete_Annotated.R` for non-interactive execution.

---

## Contributing

Enhancements, bug reports, and pull requests are welcome. Please fork the repository and submit changes via GitHub.

---

## License

This project is licensed under the **MIT License**. Include a copy of the license when distributing.
