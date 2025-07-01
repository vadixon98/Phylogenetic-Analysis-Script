############################################################
# AOS PHYLOGENY & FUNCTIONAL ANALYSIS PIPELINE
# Focused on Tomato Allene Oxide Synthase (AOS; CYP74A)
# Author: Victoria Dixon
# Date: 2025-05-01
############################################################

# === 0. Install & Load Required Packages ===

# List CRAN and Bioconductor packages needed
cran_pkgs <- c("biomartr", "ape", "phangorn", "tidyverse", "ggtree", "biomaRt")
bioc_pkgs <- c("rentrez", "msa", "Biostrings")

# Install missing CRAN packages
for(pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Install Bioconductor manager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install missing Bioconductor packages
for(pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load all libraries
library(rentrez)        # Sequence retrieval
library(biomartr)       # Alternative sequence retrieval
library(msa)            # Multiple sequence alignment
library(Biostrings)     # Sequence handling
library(ape)            # Neighbor‐Joining trees
library(phangorn)       # Maximum Likelihood phylogeny
library(ggtree)         # Tree visualization
library(tidyverse)      # Data manipulation & plotting
library(biomaRt)        # Promoter retrieval

# === 1. Sequence Retrieval ===
species_list <- c('Solanum lycopersicum', 'Solanum pimpinellifolium',
                  'Solanum pennellii', 'Solanum galapagense')
gene_name    <- 'allene oxide synthase'

# Example: retrieve CDS for reference tomato using biomartr
sl_cds <- getCodingSequence(organism = 'Solanum lycopersicum',
                            filters  = list(gene = 'CYP74A'),
                            path     = 'data')
# Write to FASTA
writeXStringSet(sl_cds, 'data/Sl_AOS_CDS.fasta')

# === 2. Promoter Retrieval ===
mart <- useEnsembl(biomart = 'plants_mart', dataset = 'slycopersicum_eg_gene')
promoters <- getSequence(id      = names(sl_cds),
                         type    = 'ensembl_gene_id',
                         seqRegion = '5utr',
                         upstream   = 2000,
                         mart    = mart)
write.fasta(as.list(promoters$sequence),
            names    = promoters$ensembl_gene_id,
            file.out = 'data/AOS_promoters_2kb.fasta')

# === 3. Multiple Sequence Alignment & Motif Scanning ===

# 3.1 Protein alignment
prot_seqs <- readAAStringSet('data/Sl_AOS_protein.fasta')
aln       <- msa(prot_seqs, method = 'ClustalW')
writeXStringSet(msaConvert(aln, type = 'bios2mds::align')[['alignment']],
                'results/AOS_protein_aln.fasta')

# 3.2 Heme‐binding motif check
pattern_heme <- 'F.GGPRC?'
heme_hits    <- vmatchPattern(pattern_heme, prot_seqs)
print(heme_hits)

# 3.3 Promoter motif scan (requires FIMO installed externally)
system('fimo --oc results/fimo_JA_elements motifs/JAElems.meme data/AOS_promoters_2kb.fasta')

# === 4. Phylogenetic Reconstruction ===

# 4.1 Read alignment and build NJ tree
dna_aln  <- read.dna('results/AOS_protein_aln.fasta', format = 'fasta')
dist_mat <- dist.dna(dna_aln, model = 'p')
NJ_tree  <- nj(dist_mat)

# 4.2 Maximum Likelihood with bootstrap
phy     <- as.phyDat(dna_aln, type = 'DNA')
fit_nj  <- pml(NJ_tree, phy)
fit_ml  <- optim.pml(fit_nj, model = 'WAG', optGamma = TRUE, optInv = TRUE)
bs      <- bootstrap.pml(fit_ml, bs = 1000, optNni = TRUE)

# 4.3 Plot trees
p_nj <- ggtree(NJ_tree) + geom_tiplab(size = 3)
p_ml <- ggtree(fit_ml$tree) +
  geom_tiplab(size = 3) +
  geom_nodelab(aes(label = round(bs, 0)), size = 2)
ggsave('results/AOS_NJ_tree.png', p_nj, width = 6, height = 8)
ggsave('results/AOS_ML_tree.png', p_ml, width = 6, height = 8)

# === 5. Functional Integration ===

# Load expression & phenotype data
expr  <- read_csv('data/AOS_expression_wound.csv')   # sample, cultivar, log2FC
pheno <- read_csv('data/AOS_resistance.csv')         # cultivar, resistance_score

# Summarize and merge
anno <- expr %>%
  group_by(cultivar) %>%
  summarize(meanFC = mean(log2FC, na.rm = TRUE)) %>%
  left_join(pheno, by = 'cultivar')

# Annotate ML tree
p_annot <- p_ml %<+% anno +
  geom_tippoint(aes(color = meanFC, size = resistance_score)) +
  scale_color_viridis_c(option = 'plasma') +
  theme_tree2()
ggsave('results/AOS_tree_annotated.png', p_annot, width = 6, height = 8)

message('Pipeline complete! Check the `results/` folder for outputs.')
