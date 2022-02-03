library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)

#-------------------------------------------------------------------------------
# Download NOMe seq data
# GCH represents accessibility scores
# HCG represents endogenous methylation scores
# https://genome.cshlp.org/content/22/12/2497.long
#-------------------------------------------------------------------------------

# HCT116
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1420nnn/GSM1420150/suppl/GSM1420150_HCT116_merge_C118BACXX_C18FYACXX_D1W2BACXX.fastq-mcf.hg19.mdups.header.clean.realign.recal.cytosine.filtered.sort.BisSNP-0.81.1.GCH.bw", 
              destfile = "path/HCT116_NOMEseq_GCH.bw")
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1420nnn/GSM1420150/suppl/GSM1420150_HCT116_merge_C118BACXX_C18FYACXX_D1W2BACXX.fastq-mcf.hg19.mdups.header.clean.realign.recal.cytosine.filtered.sort.BisSNP-0.81.1.HCG.bw", 
              destfile = "/path/HCT116_NOMEseq_HCG.bw")

# DKO1
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1420nnn/GSM1420151/suppl/GSM1420151_DKO1_merge_C118BACXX_C18G4ACXX_C1F6FACXX.hg19_rCRSchrm.fa.mdups.realign.recal.clean.cytosine.filtered.sort.BisSNP-0.81.GCH.bw", 
              destfile = "/path/DKO1_NOMEseq_GCH.bw")
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1420nnn/GSM1420151/suppl/GSM1420151_DKO1_merge_C118BACXX_C18G4ACXX_C1F6FACXX.hg19_rCRSchrm.fa.mdups.realign.recal.clean.cytosine.filtered.sort.BisSNP-0.81.HCG.bw", 
              destfile = "/path/DKO1_NOMEseq_HCG.bw")

#-------------------------------------------------------------------------------
# Read data
#-------------------------------------------------------------------------------

hct116.nome.gch = rtracklayer::import(con = "/path/to/HCT116_NOMEseq_GCH.bw")
hct116.nome.hcg = rtracklayer::import(con = "/path/to/HCT116_NOMEseq_HCG.bw")

dko1.nome.gch = rtracklayer::import(con = "/path/to/DKO1_NOMEseq_GCH.bw")
dko1.nome.hcg = rtracklayer::import(con = "/path/to/DKO1_NOMEseq_HCG.bw")

#-------------------------------------------------------------------------------
# Identify common regions between HCT116 and DKO1 so that methylated and unmethtlated regions can be defined
#-------------------------------------------------------------------------------

hct116.endo.meth <- data.frame(hct116.chr = seqnames(hct116.nome.hcg),
                               hct116.pos = ranges(hct116.nome.hcg),
                               hct116.val = hct116.nome.hcg$score)

dko1.endo.meth <- data.frame(dko1.chr = seqnames(dko1.nome.hcg),
                             dko1.pos = ranges(dko1.nome.hcg),
                             dko1.val = dko1.nome.hcg$score)

hct116.endo.meth$cmnID = stringi::stri_c(hct116.endo.meth$hct116.chr, hct116.endo.meth$hct116.pos.start, sep = "_")
dko1.endo.meth$cmnID = paste(dko1.endo.meth$dko1.chr, dko1.endo.meth$dko1.pos.start, sep = "_")
cmn.hct116.dko1_hcg  = merge(hct116.endo.meth, dko1.endo.meth, by = "cmnID")
cmn.hct116.dko1_hcg  = cmn.hct116.dko1_hcg[,2:ncol(cmn.hct116.dko1_hcg)]

# Sanity check
identical(cmn.hct116.dko1_hcg$hct116.pos.start, cmn.hct116.dko1_hcg$dko1.pos.start)
identical(as.character(cmn.hct116.dko1_hcg$hct116.chr), as.character(cmn.hct116.dko1_hcg$dko1.chr))

saveRDS(cmn.hct116.dko1_hcg, "/path/HCT116_DKO1_common_hcg.RDS")

#-------------------------------------------------------------------------------
# Extract regulatory regions
#-------------------------------------------------------------------------------

# Read reference genome
cmn.hct116.dko1_hcg <- readRDS("/path/to/HCT116_DKO1_common_hcg.RDS")
cpgrich = readRDS("/path/to/hg19_CpG_rich_promoters.RDS")
seqlevelsStyle(cpgrich) <- "UCSC"
cpgpoor <- readRDS("/path/to/hg19_CpG_poor_promoters.RDS")
seqlevelsStyle(cpgpoor) <- "UCSC"
promo <- readRDS("/path/to/hg19_promoters.RDS")
seqlevelsStyle(promo) <- "UCSC"
enh = readRDS("/path/to/hg19_enhancers.RDS")
seqlevelsStyle(enh) <- "UCSC"

# Define methylated, unMethylated, and intermediate regions
mu = cmn.hct116.dko1_hcg[cmn.hct116.dko1_hcg$hct116.val > 60 & cmn.hct116.dko1_hcg$dko1.val < 5,]
muGR = mu[,c("hct116.chr", "hct116.pos.start",  "hct116.pos.end")]
colnames(muGR) = c("chr", "start", "end")
muGR = makeGRangesFromDataFrame(muGR)

uu = cmn.hct116.dko1_hcg[cmn.hct116.dko1_hcg$hct116.val < 5 & cmn.hct116.dko1_hcg$dko1.val < 5,]
uuGR = uu[,c("hct116.chr", "hct116.pos.start",  "hct116.pos.end")]
colnames(uuGR) = c("chr", "start", "end")
uuGR = makeGRangesFromDataFrame(uuGR)

intermediate = cmn.hct116.dko1_hcg[cmn.hct116.dko1_hcg$hct116.val >= 5 & cmn.hct116.dko1_hcg$hct116.val <= 60 & cmn.hct116.dko1_hcg$dko1.val < 5,]
intermediateGR = intermediate[,c("hct116.chr", "hct116.pos.start",  "hct116.pos.end")]
colnames(intermediateGR) = c("chr", "start", "end")
intermediateGR = makeGRangesFromDataFrame(intermediateGR)

# CGI promoter
# Count how many methylated, unmethylated, and intermediate sites included in each region in reference genome
cpgrich_mu_counts <- countOverlaps(cpgrich, muGR)
cpgrich_uu_counts <- countOverlaps(cpgrich, uuGR)
cpgrich_intermediate_counts <- countOverlaps(cpgrich, intermediateGR)

if(identical(names(cpgrich_mu_counts), names(cpgrich_uu_counts))){
  # Classify into MU if MU counts are larger than 11 and proportion of MU counts are larger than 50%
  cpgrich_mu <- cpgrich[cpgrich_mu_counts > 11 & (cpgrich_mu_counts / (cpgrich_mu_counts + cpgrich_intermediate_counts + cpgrich_uu_counts)) > 0.5]
  cpgrich_uu <- cpgrich[cpgrich_uu_counts > 11 & (cpgrich_uu_counts / (cpgrich_mu_counts + cpgrich_intermediate_counts + cpgrich_uu_counts)) > 0.5]
  
  # Exclude overlapping regions between MU and UU
  cmn_cpgrich_mu <- cpgrich_mu$gene_id %in% cpgrich_uu$gene_id
  cmn_cpgrich_uu <- cpgrich_uu$gene_id %in% cpgrich_mu$gene_id
  
  cpgrich_mu <- cpgrich_mu[cmn_cpgrich_mu == F]
  cpgrich_uu <- cpgrich_uu[cmn_cpgrich_uu == F]
}

# Distal enhancer
enh_mu_counts <- countOverlaps(enh, muGR)
enh_uu_counts <- countOverlaps(enh, uuGR)
enh_intermediate_counts <- countOverlaps(enh, intermediateGR)

if(identical(names(enh_mu_counts), names(enh_uu_counts))){
  enh_mu <- enh[enh_mu_counts > 11 & (enh_mu_counts / (enh_mu_counts + enh_intermediate_counts + enh_uu_counts)) > 0.5]
  enh_uu <- enh[enh_uu_counts > 11 & (enh_uu_counts / (enh_mu_counts + enh_intermediate_counts + enh_uu_counts)) > 0.5]
  
  cmn_enh_mu <- enh_mu$genehancer_id %in% enh_uu$genehancer_id
  cmn_enh_uu <- enh_uu$genehancer_id %in% enh_mu$genehancer_id
  
  enh_mu <- enh_mu[cmn_enh_mu == F]
  enh_uu <- enh_uu[cmn_enh_uu == F]
}

# Non-CGI promoter
cpgpoor_mu_counts <- countOverlaps(cpgpoor, muGR)
cpgpoor_uu_counts <- countOverlaps(cpgpoor, uuGR)
cpgpoor_intermediate_counts <- countOverlaps(cpgpoor, intermediateGR)

if(identical(names(cpgpoor_mu_counts), names(cpgpoor_uu_counts))){
  cpgpoor_mu <- cpgpoor[cpgpoor_mu_counts > 11 & (cpgpoor_mu_counts / (cpgpoor_mu_counts + cpgpoor_intermediate_counts + cpgpoor_uu_counts)) > 0.5]
  cpgpoor_uu <- cpgpoor[cpgpoor_uu_counts > 11 & (cpgpoor_uu_counts / (cpgpoor_mu_counts + cpgpoor_intermediate_counts + cpgpoor_uu_counts)) > 0.5]
  
  cmn_cpgpoor_mu <- cpgpoor_mu$gene_id %in% cpgpoor_uu$gene_id
  cmn_cpgpoor_uu <- cpgpoor_uu$gene_id %in% cpgpoor_mu$gene_id
  
  cpgpoor_mu <- cpgpoor_mu[cmn_cpgpoor_mu == F]
  cpgpoor_uu <- cpgpoor_uu[cmn_cpgpoor_uu == F]
}

# Promoter
promo_mu_counts <- countOverlaps(promo, muGR)
promo_uu_counts <- countOverlaps(promo, uuGR)
promo_intermediate_counts <- countOverlaps(promo, intermediateGR)

if(identical(names(promo_mu_counts), names(promo_uu_counts))){
  promo_mu <- promo[promo_mu_counts > 11 & (promo_mu_counts / (promo_mu_counts + promo_intermediate_counts + promo_uu_counts)) > 0.5]
  promo_uu <- promo[promo_uu_counts > 11 & (promo_uu_counts / (promo_mu_counts + promo_intermediate_counts + promo_uu_counts)) > 0.5]
  
  cmn_promo_mu <- promo_mu$gene_id %in% promo_uu$gene_id
  cmn_promo_uu <- promo_uu$gene_id %in% promo_mu$gene_id
  
  promo_mu <- promo_mu[cmn_promo_mu == F]
  promo_uu <- promo_uu[cmn_promo_uu == F]
}

################################################################################
# Save data
################################################################################

path <- "/path/to/dir/"

if(!dir.exists(path)){
  dir.create(path, recursive = T)
}

export.bed(object = enh_mu, con = paste0(path, "enhancer_MU.bed"))
saveRDS(enh_mu, paste0(path, "enhancer_MU.RDS"))

export.bed(object = enh_uu, con = paste0(path, "enhancer_UU.bed"))
saveRDS(enh_uu, paste0(path, "enhancer_UU.RDS"))

export.bed(object = cpgrich_mu, con = paste0(path, "cpgrich_MU.bed"))
saveRDS(cpgrich_mu, paste0(path, "cpgrich_MU.RDS"))

export.bed(object = cpgrich_uu, con = paste0(path, "cpgrich_UU.bed"))
saveRDS(cpgrich_uu, paste0(path, "cpgrich_UU.RDS"))

export.bed(object = cpgpoor_mu, con = paste0(path, "cpgpoor_MU.bed"))
saveRDS(cpgpoor_mu, paste0(path, "cpgpoor_MU.RDS"))

export.bed(object = cpgpoor_uu, con = paste0(path, "cpgpoor_UU.bed"))
saveRDS(cpgpoor_uu, paste0(path, "cpgpoor_UU.RDS"))

export.bed(object = promo_mu, con = paste0(path, "promo_MU.bed"))
saveRDS(promo_mu, paste0(path, "promo_MU.RDS"))

export.bed(object = promo_uu, con = paste0(path, "promo_UU.bed"))
saveRDS(promo_uu, paste0(path, "promo_UU.RDS"))

export.bed(object = cpgrich, con = paste0(path, "cpgrich.bed"))
saveRDS(cpgrich, paste0(path, "cpgrich.RDS"))

export.bed(object = cpgpoor, con = paste0(path, "cpgpoor.bed"))
saveRDS(cpgpoor, paste0(path, "cpgpoor.RDS"))
