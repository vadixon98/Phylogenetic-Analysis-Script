############################################################
# AOS PHYLOGENY & FUNCTIONAL ANALYSIS PIPELINE
# Focused on Tomato Allene Oxide Synthase (AOS; CYP74A)
# Author: Victoria Dixon
# Date: 2025-05-01
############################################################

# === 0. Install & Load Required Packages ===
# CRAN packages needed
cran_pkgs <- c("biomartr", "ape", "phangorn", "tidyverse", "seqinr")
# Bioconductor packages needed
bioc_pkgs <- c("rentrez", "msa", "Biostrings", "ggtree", "biomaRt")

# Install missing CRAN packages
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install missing Bioconductor packages
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load all libraries
invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

# === 1. Sequence Retrieval ===
message('Downloading CDS and proteome for Solanum lycopersicum...')
# Fetch and save CDS
cds_file <- getCDS(organism = 'Solanum lycopersicum', path = 'data')
all_cds  <- Biostrings::readDNAStringSet(cds_file)
# Fetch and save proteome
prot_file <- getProteome(organism = 'Solanum lycopersicum', path = 'data')
all_prot <- Biostrings::readAAStringSet(prot_file)

# Subset for AOS (CYP74A) entries by header pattern
aos_cds  <- all_cds[grep('CYP74A', names(all_cds))]
aos_prot <- all_prot[grep('CYP74A', names(all_prot))]

# Write filtered FASTA files
Biostrings::writeXStringSet(aos_cds,  'data/Sl_AOS_CDS.fasta')
Biostrings::writeXStringSet(aos_prot, 'data/Sl_AOS_protein.fasta')

# === 2. Promoter Retrieval ===
message('Automated promoter retrieval via Ensembl Plants BioMart...')
# List marts and datasets for debugging
print(biomaRt::listMarts())
plants_mart <- biomaRt::useMart(
  biomart = 'plants_mart',
  host    = 'https://plants.ensembl.org'
)
print(head(biomaRt::listDatasets(plants_mart)))
# Connect to the 'slycopersicum_eg_gene' dataset
mart2 <- biomaRt::useMart(
  biomart = 'plants_mart',
  dataset = 'slycopersicum_eg_gene',
  host    = 'https://plants.ensembl.org'
)
# Prepare transcript and gene IDs
trans_ids <- sub("\\..*", "", names(aos_cds))
# Attempt mapping transcripts to genes
map_df <- biomaRt::getBM(
  attributes = c('ensembl_transcript_id','ensembl_gene_id'),
  filters    = 'ensembl_transcript_id',
  values     = trans_ids,
  mart       = mart2
)
gene_ids <- unique(map_df$ensembl_gene_id)

promoter_fasta <- 'data/AOS_promoters_2kb.fasta'
promoter_success <- FALSE

# Primary method: getBM upstream_gene_flank
if(length(gene_ids)>0) {
  try({
    prom_df <- biomaRt::getBM(
      attributes = c('ensembl_gene_id','upstream_gene_flank'),
      filters    = 'ensembl_gene_id',
      values     = gene_ids,
      mart       = mart2
    )
    if(nrow(prom_df)>0) {
      seqinr::write.fasta(
        as.list(prom_df$upstream_gene_flank),
        names = prom_df$ensembl_gene_id,
        file.out = promoter_fasta
      )
      promoter_success <- TRUE
      message('Promoters retrieved via getBM upstream_gene_flank')
    }
  }, silent=TRUE)
}
# Fallback method: biomartr::getSequence type='upstream'
if(!promoter_success && length(gene_ids)>0) {
  message('Falling back to biomartr::getSequence for upstream sequences...')
  try({
    seqs2 <- biomartr::getSequence(
      db        = 'plants_mart',
      organism  = 'Solanum lycopersicum',
      type      = 'upstream',
      id        = gene_ids,
      up        = 2000,
      download  = FALSE
    )
    seqinr::write.fasta(
      as.list(seqs2$sequence),
      names = seqs2$id,
      file.out = promoter_fasta
    )
    promoter_success <- TRUE
    message('Promoters retrieved via biomartr::getSequence')
  }, silent=TRUE)
}
if(!promoter_success) {
  stop('Automated promoter retrieval failed: please check dataset names or network')
}

# === 3. Multiple Sequence Alignment & Motif Scanning === Multiple Sequence Alignment & Motif Scanning === Multiple Sequence Alignment & Motif Scanning === Multiple Sequence Alignment & Motif Scanning ===
# 3.1 Protein alignment
prot_seqs <- Biostrings::readAAStringSet('data/Sl_AOS_protein.fasta')
aln       <- msa::msa(prot_seqs, method = 'ClustalW')
# Export alignment
msa::msaPrettyPrint(aln, output = 'pdf', file = 'results/AOS_protein_aln.pdf')

# 3.2 Heme-binding motif check
pattern_heme <- 'F.GGPRC?'
heme_hits    <- Biostrings::vmatchPattern(pattern_heme, prot_seqs)
print(heme_hits)

# 3.3 Promoter motif scan (requires FIMO externally)
system('fimo --oc results/fimo_JA_elements motifs/JAElems.meme data/AOS_promoters_2kb.fasta')

# === 4. Phylogenetic Reconstruction ===
# Read amino-acid alignment for tree building (convert to DNA alignment if needed)
dna_aln  <- ape::read.dna('results/AOS_protein_aln.fasta', format = 'fasta')
dist_mat <- ape::dist.dna(dna_aln, model = 'p')
NJ_tree  <- ape::nj(dist_mat)

# Maximum Likelihood reconstruction
phy    <- phangorn::as.phyDat(dna_aln, type = 'DNA')
fit_nj <- phangorn::pml(NJ_tree, phy)
fit_ml <- phangorn::optim.pml(fit_nj, model = 'WAG', optGamma = TRUE, optInv = TRUE)
bs     <- phangorn::bootstrap.pml(fit_ml, bs = 1000, optNni = TRUE)

# Plot and save trees
gg_nj <- ggtree::ggtree(NJ_tree) + ggtree::geom_tiplab(size = 3)
ggsave('results/AOS_NJ_tree.png', gg_nj, width = 6, height = 8)

gg_ml <- ggtree::ggtree(fit_ml$tree) +
  ggtree::geom_tiplab(size = 3) +
  ggtree::geom_nodelab(aes(label = round(bs, 0)), size = 2)
ggsave('results/AOS_ML_tree.png', gg_ml, width = 6, height = 8)

# === 5. Functional Integration ===
expr  <- readr::read_csv('data/AOS_expression_wound.csv')   # sample, cultivar, log2FC
pheno <- readr::read_csv('data/AOS_resistance.csv')         # cultivar, resistance_score

df_anno <- expr %>%
  dplyr::group_by(cultivar) %>%
  dplyr::summarize(meanFC = mean(log2FC, na.rm = TRUE)) %>%
  dplyr::left_join(pheno, by = 'cultivar')

annotated_tree <- gg_ml %<+% df_anno +
  ggtree::geom_tippoint(aes(color = meanFC, size = resistance_score)) +
  ggplot2::scale_color_viridis_c(option = 'plasma') +
  ggtree::theme_tree2()

ggsave('results/AOS_tree_annotated.png', annotated_tree, width = 6, height = 8)

message('Pipeline complete! Check the `data/` and `results/` folders for outputs.')
