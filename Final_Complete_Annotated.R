############################################################
# AOS PHYLOGENY & FUNCTIONAL ANALYSIS PIPELINE
# Script: Final_Complete_Annotated.R
# Description: This pipeline retrieves Allene Oxide Synthase (AOS; CYP74A) sequences,
#              extracts promoter regions, performs multiple sequence alignment,
#              scans for motifs, and builds phylogenetic trees with bootstrap support.
# Author: Victoria A Dixon
# Date: 2025-05-02 
############################################################

# === Setup Directories ===
# Create 'data' and 'results' directories if they do not exist to organize outputs
# The 'data' directory will store downloaded sequences, genomes, and intermediate files
# The 'results' directory will contain alignments, trees, and analysis outputs
if (!dir.exists('data')) dir.create('data', recursive = TRUE)
if (!dir.exists('results')) dir.create('results', recursive = TRUE)

# === 0. Install & Load Required Packages ===
# Define CRAN and Bioconductor package lists
# CRAN packages: biomartr (genome data), ape/phangorn (phylogenetics), 
#                tidyverse (data manipulation), seqinr (sequence utilities)
cran_pkgs <- c("biomartr", "ape", "phangorn", "tidyverse", "seqinr")
# Bioconductor packages: rentrez (NCBI access), msa (multiple alignment),
#                        Biostrings (sequence handling), ggtree (tree visualization),
#                        biomaRt/rtracklayer/IRanges (genome annotation)
bioc_pkgs <- c("rentrez", "msa", "Biostrings", "ggtree", "biomaRt", "rtracklayer", "IRanges")

# Install missing CRAN packages from specified repository
# Checks each package and only installs if not already available
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Ensure BiocManager is available to install Bioconductor packages
# BiocManager is required to install packages from Bioconductor repositories
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install missing Bioconductor packages using BiocManager
# Bioconductor packages require a different installation method than CRAN
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load all required libraries into the session
# invisible() suppresses output from lapply; character.only = TRUE allows string package names
invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

# === 1. Sequence Retrieval ===
# Download coding sequences (CDS) and protein sequences for tomato
# This section retrieves all coding sequences and proteins from NCBI/Ensembl
# and filters for Allene Oxide Synthase (AOS) genes
message('Downloading CDS and proteome for Solanum lycopersicum...')
# Download all coding sequences (CDS) for tomato and save to data directory
cds_file  <- getCDS(organism = 'Solanum lycopersicum', path = 'data')
# Read CDS sequences into a DNAStringSet object for manipulation
all_cds   <- Biostrings::readDNAStringSet(cds_file)
# Download all protein sequences for tomato and save to data directory
prot_file <- getProteome(organism = 'Solanum lycopersicum', path = 'data')
# Read protein sequences into an AAStringSet object for manipulation
all_prot  <- Biostrings::readAAStringSet(prot_file)

# Define pattern to identify AOS sequences in headers (case-insensitive)
# Pattern matches: CYP74A (gene family), "allene oxide synthase" (full name), or "AOS" (abbreviation)
pattern   <- 'CYP74A|allene oxide synthase|AOS'
# Subset sequences matching the AOS pattern in sequence names/headers
# grep() searches for pattern matches; ignore.case = TRUE makes search case-insensitive
aos_cds   <- all_cds[grep(pattern, names(all_cds), ignore.case = TRUE)]
aos_prot  <- all_prot[grep(pattern, names(all_prot), ignore.case = TRUE)]

# Check that sequences were found, else stop with informative error
# This prevents downstream errors if no AOS sequences are present in the dataset
if (length(aos_cds) == 0) stop('No AOS CDS sequences found; check header patterns')
if (length(aos_prot) == 0) stop('No AOS protein sequences found; check header patterns')

# Write selected AOS sequences to FASTA files for downstream analysis
# These files will be used for alignment and phylogenetic analysis
Biostrings::writeXStringSet(aos_cds,  'data/Sl_AOS_CDS.fasta')
Biostrings::writeXStringSet(aos_prot, 'data/Sl_AOS_protein.fasta')

# === 2. Promoter Retrieval ===
# Extract ~2kb upstream promoter sequences for AOS genes
# Promoters are important for understanding gene regulation and contain cis-regulatory elements
message('Retrieving ~2kb promoter sequences using GFF and genome FASTA...')
# Download GFF3 annotation file (contains gene coordinates and annotations)
gff3_file <- biomartr::getGFF(organism = 'Solanum lycopersicum', path = 'data')
# Download complete genome FASTA file (contains all chromosome sequences)
genome_fa <- biomartr::getGenome(organism = 'Solanum lycopersicum', path = 'data')
# Import GFF3 file as a GRanges object for coordinate-based operations
gff       <- rtracklayer::import.gff3(gff3_file)
# Read genome sequences into a DNAStringSet for sequence extraction
genome    <- Biostrings::readDNAStringSet(genome_fa)

# Extract gene features matching the AOS pattern in annotations
# mcols() extracts metadata columns from GRanges object
attrs     <- mcols(gff)$attributes
# Filter GFF to only AOS genes: type == 'gene' and pattern match in attributes
# grepl() searches for pattern in the attributes string (gene descriptions)
aos_genes <- gff[gff$type == 'gene' & grepl('allene oxide synthase|CYP74A|AOS', attrs, ignore.case = TRUE)]
# Error check: ensure at least one AOS gene was found
if (length(aos_genes) == 0) stop('No AOS gene found; check GFF attributes')

# Function to extract ~2kb upstream promoter sequence for each gene
# lapply() applies the function to each AOS gene found
promoter_seqs <- lapply(aos_genes, function(g) {
  # Extract chromosome name and strand orientation for the gene
  chr    <- as.character(seqnames(g))
  strand <- as.character(strand(g))
  # Determine start and end positions based on gene orientation
  # For + strand: promoter is upstream (before start), so subtract 2000bp
  if (strand == '+') {
    st <- start(g) - 2000; ed <- start(g) - 1
  } else {
    # For - strand: promoter is downstream (after end), so add 2000bp
    # Note: coordinates will be reversed later via reverse-complement
    st <- end(g) + 1; ed <- end(g) + 2000
  }
  # Ensure coordinates do not fall below 1 (chromosome start)
  st <- max(st, 1)
  # Extract sequence from genome using chromosome name and coordinates
  seq <- Biostrings::getSeq(genome, names = chr, start = st, end = ed)
  # Convert to DNAStringSet format
  seq_set <- Biostrings::DNAStringSet(as.character(seq))
  # Reverse-complement if gene is on the minus strand
  # This ensures all promoters are in 5' to 3' orientation relative to transcription start
  if (strand == '-') seq_set <- Biostrings::reverseComplement(seq_set)
  # Name sequence by gene ID or Name attribute (whichever is available)
  nm <- if (!is.null(g$Name) && nzchar(g$Name)) g$Name else g$ID
  names(seq_set) <- nm
  seq_set
})
# Combine promoter sequences from all genes into a single DNAStringSet
# Extract character strings from each DNAStringSet in the list
promoter_chars <- sapply(promoter_seqs, function(x) as.character(x)[1])
# Create a unified DNAStringSet from all promoter sequences
promoters_set  <- Biostrings::DNAStringSet(promoter_chars)
# Preserve gene names from the original list
names(promoters_set) <- names(promoter_seqs)
# Save promoters to FASTA for motif analysis (e.g., FIMO scanning)
Biostrings::writeXStringSet(promoters_set, 'data/AOS_promoters_2kb.fasta')

# === 3. Multiple Sequence Alignment & Motif Scanning ===
# Align AOS protein sequences to identify conserved regions and evolutionary relationships
# ClustalW is a progressive alignment algorithm suitable for protein sequences
prot_seqs <- Biostrings::readAAStringSet('data/Sl_AOS_protein.fasta')
# Perform multiple sequence alignment using ClustalW algorithm
# ClustalW uses progressive alignment: aligns most similar sequences first, then adds others
aln       <- msa::msa(prot_seqs, method = 'ClustalW')
# Convert alignment object to AAStringSet format for saving and further analysis
msa_aa    <- as(aln, "AAStringSet")
# Save alignment for record and downstream phylogenetic analysis
Biostrings::writeXStringSet(msa_aa, 'results/AOS_protein_aln.fasta')

# Example motif scan: heme-binding motif F-GGPRC in protein sequences
# AOS is a cytochrome P450 enzyme; heme-binding domain is essential for function
# Pattern uses '.' as a wildcard to match any single amino acid
pattern_heme <- 'F.GGPRC'
# Search for motif pattern in all protein sequences
# vmatchPattern() returns match positions for each sequence
heme_hits    <- Biostrings::vmatchPattern(pattern_heme, prot_seqs)
# Display motif hit positions for user review
# This helps verify that AOS sequences contain expected functional domains
print(heme_hits)

# Run FIMO from the MEME suite to scan promoter sequences for jasmonate-responsive elements
# FIMO (Find Individual Motif Occurrences) scans sequences against a motif database
# JAElems.meme should contain jasmonate-responsive cis-regulatory elements (e.g., G-box, GCC-box)
# Note: Requires FIMO to be installed separately (part of MEME suite)
# --oc specifies output directory; results will contain motif matches with p-values
system('fimo --oc results/fimo_JA_elements motifs/JAElems.meme data/AOS_promoters_2kb.fasta')

# === 4. Phylogenetic Reconstruction ===
# Build phylogenetic trees to infer evolutionary relationships among AOS sequences
# Uses Maximum Likelihood (ML) method with bootstrap support for node confidence
message('Building phylogeny from protein alignment...')
# Convert MSA to AAbin format (ape package format for amino acid sequences)
# AAbin is required for distance calculations and tree building in ape/phangorn
alignment_aabin <- msa::msaConvert(aln, type = "ape::AAbin")
# Convert to phyDat format (phangorn's internal format for phylogenetic data)
# type = 'AA' specifies amino acid data (vs 'DNA' for nucleotide)
phydat          <- phangorn::phyDat(alignment_aabin, type = 'AA')

# Calculate maximum likelihood tree with JTT model and gamma distribution
# JTT (Jones-Taylor-Thornton) is a standard substitution model for protein evolution
# Gamma distribution accounts for rate variation across sites
# dist.ml() calculates pairwise distances using maximum likelihood
dist_mat <- phangorn::dist.ml(phydat)
# Build Neighbor-Joining (NJ) tree as a starting point for ML optimization
# NJ is fast and provides a reasonable initial tree topology
NJ_tree  <- ape::nj(dist_mat)
# Fit ML model to NJ tree to get initial likelihood score
# pml() creates a pml object with tree, data, and model parameters
fit_nj   <- phangorn::pml(NJ_tree, phydat)
# Optimize ML tree: adjust branch lengths, estimate gamma shape parameter, and proportion of invariable sites
# optGamma = TRUE estimates gamma distribution shape parameter
# optInv = TRUE estimates proportion of invariable sites
# This finds the tree topology and parameters that maximize likelihood
fit_ml   <- phangorn::optim.pml(fit_nj, model = 'JTT', optGamma = TRUE, optInv = TRUE)

# Bootstrap analysis to assess node support (1000 replicates)
# Bootstrap resamples alignment columns and rebuilds trees to test node stability
# bs = 1000 means 1000 bootstrap replicates (higher = more reliable but slower)
# optNni = TRUE optimizes tree topology during each bootstrap (nearest-neighbor interchanges)
trees_bs <- phangorn::bootstrap.pml(fit_ml, bs = 1000, optNni = TRUE)
# Transfer bootstrap support values to the ML tree nodes
# plotBS() adds bootstrap percentages as node labels (quiet = TRUE suppresses printing)
bs_tree   <- phangorn::plotBS(fit_ml$tree, trees_bs, quiet = TRUE)

# Plot and save Neighbor-Joining tree with tip labels
# ggtree provides ggplot2-based tree visualization
# geom_tiplab() adds sequence names at tree tips
gg_nj <- ggtree::ggtree(NJ_tree) + ggtree::geom_tiplab(size = 3)
# Save tree as PNG image (width/height in inches)
ggsave('results/AOS_NJ_tree.png', gg_nj, width = 6, height = 8)

# Plot and save ML tree with bootstrap support labels on nodes
# bs_tree contains bootstrap support values as node labels
# geom_nodelab() adds bootstrap percentages to internal nodes
# hjust = -0.3 positions labels slightly to the left of nodes
p_ml <- ggtree::ggtree(bs_tree) + 
  ggtree::geom_tiplab(size = 3) +
  ggtree::geom_nodelab(aes(label = label), hjust = -0.3, size = 2)
# Save ML tree with bootstrap support values
ggsave('results/AOS_ML_tree.png', p_ml, width = 6, height = 8)

# === 5. Pipeline Completion ===
# All analysis steps completed successfully
# Output files include:
#   - data/: AOS sequences (CDS, protein), promoters, genome files
#   - results/: Alignments, phylogenetic trees (NJ and ML), FIMO motif scan results
message('Pipeline complete! Results are saved in the results directory.')
