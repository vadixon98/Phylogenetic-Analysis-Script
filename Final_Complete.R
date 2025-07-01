############################################################
# AOS PHYLOGENY & FUNCTIONAL ANALYSIS PIPELINE
# Focused on Tomato Allene Oxide Synthase (AOS; CYP74A)
# Author: Victoria Dixon
# Date: 2025-05-01
############################################################

# === Setup Directories ===
if (!dir.exists('data')) dir.create('data', recursive=TRUE)
if (!dir.exists('results')) dir.create('results', recursive=TRUE)

# === 0. Install & Load Required Packages ===
cran_pkgs <- c("biomartr", "ape", "phangorn", "tidyverse", "seqinr")
bioc_pkgs <- c("rentrez", "msa", "Biostrings", "ggtree", "biomaRt", "rtracklayer", "IRanges")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg, repos="https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg)
}

invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only=TRUE))

# === 1. Sequence Retrieval ===
message('Downloading CDS and proteome for Solanum lycopersicum...')
cds_file  <- getCDS(organism='Solanum lycopersicum', path='data')
all_cds   <- Biostrings::readDNAStringSet(cds_file)
prot_file <- getProteome(organism='Solanum lycopersicum', path='data')
all_prot  <- Biostrings::readAAStringSet(prot_file)

pattern   <- 'CYP74A|allene oxide synthase|AOS'
aos_cds   <- all_cds[grep(pattern, names(all_cds), ignore.case=TRUE)]
aos_prot  <- all_prot[grep(pattern, names(all_prot), ignore.case=TRUE)]
if (length(aos_cds)==0) stop('No AOS CDS sequences found; check header patterns')
if (length(aos_prot)==0) stop('No AOS protein sequences found; check header patterns')

Biostrings::writeXStringSet(aos_cds,  'data/Sl_AOS_CDS.fasta')
Biostrings::writeXStringSet(aos_prot, 'data/Sl_AOS_protein.fasta')

# === 2. Promoter Retrieval ===
message('Retrieving ~2kb promoter sequences using GFF and genome FASTA...')
gff3_file <- biomartr::getGFF(organism='Solanum lycopersicum', path='data')
genome_fa <- biomartr::getGenome(organism='Solanum lycopersicum', path='data')
gff       <- rtracklayer::import.gff3(gff3_file)
genome    <- Biostrings::readDNAStringSet(genome_fa)
attrs     <- mcols(gff)$attributes
aos_genes <- gff[gff$type=='gene' & grepl('allene oxide synthase|CYP74A|AOS', attrs, ignore.case=TRUE)]
if (length(aos_genes)==0) stop('No AOS gene found; check GFF attributes')

promoter_seqs <- lapply(aos_genes, function(g) {
  chr    <- as.character(seqnames(g))
  strand <- as.character(strand(g))
  if (strand=='+') { st <- start(g)-2000; ed <- start(g)-1 } else { st <- end(g)+1; ed <- end(g)+2000 }
  st <- max(st, 1)
  seq <- Biostrings::getSeq(genome, names=chr, start=st, end=ed)
  seq_set <- Biostrings::DNAStringSet(as.character(seq))
  if (strand=='-') seq_set <- Biostrings::reverseComplement(seq_set)
  nm <- if (!is.null(g$Name) && nzchar(g$Name)) g$Name else g$ID
  names(seq_set) <- nm
  seq_set
})
promoter_chars <- sapply(promoter_seqs, function(x) as.character(x)[1])
promoters_set  <- Biostrings::DNAStringSet(promoter_chars)
names(promoters_set) <- names(promoter_seqs)
Biostrings::writeXStringSet(promoters_set, 'data/AOS_promoters_2kb.fasta')

# === 3. Multiple Sequence Alignment & Motif Scanning ===
prot_seqs <- Biostrings::readAAStringSet('data/Sl_AOS_protein.fasta')
aln       <- msa::msa(prot_seqs, method='ClustalW')
msa_aa    <- as(aln, "AAStringSet")
Biostrings::writeXStringSet(msa_aa, 'results/AOS_protein_aln.fasta')

pattern_heme <- 'F.GGPRC'
heme_hits    <- Biostrings::vmatchPattern(pattern_heme, prot_seqs)
print(heme_hits)

system('fimo --oc results/fimo_JA_elements motifs/JAElems.meme data/AOS_promoters_2kb.fasta')

# === 4. Phylogenetic Reconstruction ===
message('Building phylogeny from protein alignment...')
# Convert alignment to AAbin and then to phyDat
alignment_aabin <- msa::msaConvert(aln, type = "ape::AAbin")
phydat          <- phangorn::phyDat(alignment_aabin, type='AA')
# Compute ML tree with bootstrap support
dist_mat <- phangorn::dist.ml(phydat)
NJ_tree  <- ape::nj(dist_mat)
fit_nj   <- phangorn::pml(NJ_tree, phydat)
fit_ml   <- phangorn::optim.pml(fit_nj, model='JTT', optGamma=TRUE, optInv=TRUE)
# Perform bootstrap replicates
trees_bs <- phangorn::bootstrap.pml(fit_ml, bs=1000, optNni=TRUE)
# Annotate ML tree with bootstrap support using plotBS
bs_tree <- phangorn::plotBS(fit_ml$tree, trees_bs, quiet = TRUE)

# Plot Neighbor-Joining tree
gg_nj <- ggtree::ggtree(NJ_tree) + ggtree::geom_tiplab(size=3)
ggsave('results/AOS_NJ_tree.png', gg_nj, width=6, height=8)

# Plot ML tree with bootstrap support labels
p_ml <- ggtree::ggtree(bs_tree) + 
  ggtree::geom_tiplab(size=3) +
  ggtree::geom_nodelab(aes(label=label), hjust=-0.3, size=2)
ggsave('results/AOS_ML_tree.png', p_ml, width=6, height=8)

# === 5. Complete ===
message('Pipeline complete!')
