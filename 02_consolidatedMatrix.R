library(GenomicRanges)
library(recount)
library(rtracklayer)

path = "/path/to/dir/" # Access to data (e.g., BED files for MU, UU, or methylation regions or expression data)
path2 = "/path/to/dir/"

#-------------------------------------------------------------------------------
# Read regulatory MU and UU regions
#-------------------------------------------------------------------------------

cpgrich_MU = readRDS(paste0(path, "cpgrich_MU.RDS"))
cpgrich_UU = readRDS(paste0(path, "cpgrich_UU.RDS"))
# enhancer_MU = readRDS(paste0(path, "enhancer_MU.RDS"))
# enhancer_UU = readRDS(paste0(path, "enhancer_UU.RDS"))
cpgpoor_MU = readRDS(paste0(path, "cpgpoor_MU.RDS"))
cpgpoor_UU = readRDS(paste0(path, "cpgpoor_UU.RDS"))
# colnames(mcols(enhancer_MU))[6] = "gene_name"
# colnames(mcols(enhancer_UU))[6] = "gene_name"

#-------------------------------------------------------------------------------
# Read HCT116 histone marks and chromatin accessibility (ATAC)
#-------------------------------------------------------------------------------

cleanData = function(file){
  tab = data.table::fread(file)
  colnames(tab) = gsub("'", "", colnames(tab))
  colnames(tab) = gsub("#", "", colnames(tab))
  tab = data.frame(tab, stringsAsFactors = F)
  tab = makeGRangesFromDataFrame(tab, keep.extra.columns = T)
  return(tab)
}

# Access to summarized score matrices by deepTools
HCT116_cpgrich_MU = cleanData(file = paste0(path2, "multiBigwigSummary/HCT116/tabular/cpgrich_MU_summarized_scores.tab"))
HCT116_cpgrich_UU = cleanData(file = paste0(path2, "multiBigwigSummary/HCT116/tabular/cpgrich_UU_summarized_scores.tab"))
# HCT116_enhancer_MU = cleanData(file = paste0(path2, "multiBigwigSummary/HCT116/tabular/enhancer_MU_summarized_scores.tab"))
# HCT116_enhancer_UU = cleanData(file = paste0(path2, "multiBigwigSummary/HCT116/tabular/enhancer_UU_summarized_scores.tab"))
HCT116_cpgpoor_MU = cleanData(file = paste0(path2, "multiBigwigSummary/HCT116/tabular/cpgpoor_MU_summarized_scores.tab"))
HCT116_cpgpoor_UU = cleanData(file = paste0(path2, "multiBigwigSummary/HCT116/tabular/cpgpoor_UU_summarized_scores.tab"))
rm(cleanData)

#-------------------------------------------------------------------------------
# Read DKO1 histone marks and chromatin accessibility (ATAC)
#-------------------------------------------------------------------------------

cleanData = function(file){
  tab = data.table::fread(file)
  colnames(tab) = gsub("'", "", colnames(tab))
  colnames(tab) = gsub("#", "", colnames(tab))
  tab = data.frame(tab, stringsAsFactors = F)

  # Removing the bigwig summarized expression values. 
  # We will use the recount computed TPM values (see below)
  tab = tab[, -grep("Expression", colnames(tab))]
  
  #for comparing Bigwig expression and Recount expression
  #colnames(tab) <- gsub(".*Expression_Rep1.*", "BW_expression_Rep1", colnames(tab))
  #colnames(tab) <- gsub(".*Expression_Rep2.*", "BW_expression_Rep2", colnames(tab))
  
  tab = makeGRangesFromDataFrame(data.frame(tab, stringsAsFactors = F), keep.extra.columns = T)
  return(tab)
}

DKO1_cpgrich_MU = cleanData(file = paste0(path2, "multiBigwigSummary/DKO1/tabular/cpgrich_MU_summarized_scores.tab"))
DKO1_cpgrich_UU = cleanData(file = paste0(path2, "multiBigwigSummary/DKO1/tabular/cpgrich_UU_summarized_scores.tab"))
# DKO1_enhancer_MU = cleanData(file = paste0(path2, "multiBigwigSummary/DKO1/tabular/enhancer_MU_summarized_scores.tab"))
# DKO1_enhancer_UU = cleanData(file = paste0(path2, "multiBigwigSummary/DKO1/tabular/enhancer_UU_summarized_scores.tab"))
DKO1_cpgpoor_MU = cleanData(file = paste0(path2, "multiBigwigSummary/DKO1/tabular/cpgpoor_MU_summarized_scores.tab"))
DKO1_cpgpoor_UU = cleanData(file = paste0(path2, "multiBigwigSummary/DKO1/tabular/cpgpoor_UU_summarized_scores.tab"))
rm(cleanData)

#-------------------------------------------------------------------------------
# Sanity check if HCT116 and DKO1 MU/UU share the same regions
#-------------------------------------------------------------------------------

identical(length(ranges(HCT116_cpgrich_MU)), length(ranges(DKO1_cpgrich_MU)))
identical(length(ranges(HCT116_cpgrich_UU)), length(ranges(DKO1_cpgrich_UU)))
# identical(length(ranges(HCT116_enhancer_MU)), length(ranges(DKO1_enhancer_MU)))
# identical(length(ranges(HCT116_enhancer_UU)), length(ranges(DKO1_enhancer_UU)))
identical(length(ranges(HCT116_cpgpoor_MU)), length(ranges(DKO1_cpgpoor_MU)))
identical(length(ranges(HCT116_cpgpoor_UU)), length(ranges(DKO1_cpgpoor_UU)))

#-------------------------------------------------------------------------------
# Read HCT116 methylation
#-------------------------------------------------------------------------------

HCT116_meth = data.table::fread(paste0(path, "HCT116/HCT116_GSM1465024_methylation.bed.gz"), skip = 1)
cor(HCT116_meth$V5, HCT116_meth$V7)
#[1] 0.9999998
HCT116_meth = data.frame(HCT116_meth[,c(1,2,3,7)], stringsAsFactors = F)
colnames(HCT116_meth) = c("chr", "start", "end", "methval")
HCT116_meth = makeGRangesFromDataFrame(HCT116_meth, keep.extra.columns = T)

#-------------------------------------------------------------------------------
# Read DKO1 methylation
#-------------------------------------------------------------------------------

DKO1_meth = data.table::fread(paste0(path, "DKO1/DKO1_GSE60106_methylation.bed.gz"), skip = 1)
cor(DKO1_meth$V5, DKO1_meth$V7)
#[1] 0.9999997
DKO1_meth = data.frame(DKO1_meth[,c(1,2,3,7)], stringsAsFactors = F)
colnames(DKO1_meth) = c("chr", "start", "end", "methval")
DKO1_meth = makeGRangesFromDataFrame(DKO1_meth, keep.extra.columns = T)

#-------------------------------------------------------------------------------
# Read HCT116 RNAseq
#-------------------------------------------------------------------------------

exp1 = read.table(paste0(path, "HCT116/HCT116_GSE52429_FPKM_repPF.txt.gz"), header=T, stringsAsFactors = F)
exp1$gene_id <- gsub("\\..*", "", exp1$gene_id)
exp1 = data.frame(Gene_id = exp1$gene_id ,Gene = exp1$gene_short_name, Expression_Rep1 = exp1$FPKM, stringsAsFactors = F)

exp2 = read.table(paste0(path, "HCT116/HCT116_GSE52429_FPKM_repPJ.txt.gz"), header=T, stringsAsFactors = F)
exp2$gene_id <- gsub("\\..*", "", exp2$gene_id)
exp2 = data.frame(Gene_id = exp2$gene_id ,Gene = exp2$gene_short_name, Expression_Rep2 = exp2$FPKM, stringsAsFactors = F)

GROI <- list(cpgrich_MU, cpgrich_UU, cpgpoor_MU, cpgpoor_UU) # enhancer_MU, enhancer_UU
HCT116_exp_list <- list()
HCT116_exp_corr_list <- list()
for(i in 1:length(GROI)){
  exp_a = exp1[exp1$Gene_id %in% GROI[[i]]$gene_id,]
  exp_b = exp2[exp2$Gene_id %in% GROI[[i]]$gene_id,]
  exp = merge(exp_a, exp_b, by = "Gene_id")

  # Take advantage of gene names since gene ID do not match (adjust any code if needed down below)
  # exp_a = exp1[exp1$Gene %in% enhancer$gene_name,]
  # exp_a = aggregate(Expression_Rep1 ~ Gene, data = exp_a, FUN = max)
  # exp_b = exp2[exp2$Gene %in% enhancer$gene_name,]
  # exp_b = aggregate(Expression_Rep2 ~ Gene, data = exp_b, FUN = max)
  # hct116_exp = merge(exp_a, exp_b)

  HCT116_exp_list[[i]] <- exp
  corr <- cor(exp$Expression_Rep1, exp$Expression_Rep2, method = "spearman")
  HCT116_exp_corr_list[[i]] <- corr
}

rm(exp1, exp2, exp_a, exp_b, exp, GROI)
names(HCT116_exp_list) <- c("cpgrich_MU", "cpgrich_UU", "cpgpoor_MU", "cpgpoor_UU") # "enhancer_MU", "enhancer_UU"
names(HCT116_exp_corr_list) <- c("cpgrich_MU", "cpgrich_UU", "cpgpoor_MU", "cpgpoor_UU") # "enhancer_MU", "enhancer_UU"

#-------------------------------------------------------------------------------
# DKO1 RNAseq
#-------------------------------------------------------------------------------

project_info = abstract_search("DKO1")
download_study(project_info$project)
load(file.path(project_info$project, "rse_gene.Rdata"))
unlink(project_info$project, recursive = TRUE)
exp = getTPM(rse_gene)

# remove those whose symbol is NA
if(identical(rowData(rse_gene)$gene_id, rownames(exp))){
  g = rowData(rse_gene)$symbol
  rmv = unique(c(which(unlist(sapply(g, is.na))), which(sapply(g, length) > 1)))
  exp = exp[- rmv,]
  g = rowData(rse_gene)
  g = g[- rmv,]
  
  genanno = data.frame(gene_id = unlist(g$gene_id),
                       bp_length = unlist(g$bp_length),
                       symbol = unlist(g$symbol), 
                       stringsAsFactors = F)
  rm(g, rmv)
}

if(identical(genanno$gene_id, rownames(exp))){
  genanno$gene_id <- gsub("\\..*", "", genanno$gene_id)
  exp1 = data.frame(Gene_id = genanno$gene_id, Gene = genanno$symbol, Expression_Rep1 = exp[,1], stringsAsFactors = F)
  exp2 = data.frame(Gene_id = genanno$gene_id, Gene = genanno$symbol, Expression_Rep2 = exp[,2], stringsAsFactors = F)
  
  GROI <- list(cpgrich_MU, cpgrich_UU, cpgpoor_MU, cpgpoor_UU) # enhancer_MU, enhancer_UU
  DKO1_exp_list <- list()
  DKO1_exp_corr_list <- list()
  for(i in 1:length(GROI)){
    exp_a = exp1[exp1$Gene_id %in% GROI[[i]]$gene_id,]
    exp_b = exp2[exp2$Gene_id %in% GROI[[i]]$gene_id,]
    exp = merge(exp_a, exp_b, "Gene_id")

    # Take advantage of gene names since gene ID do not match
    # exp_a = exp1[exp1$Gene %in% enhancer$gene_name,]
    # exp_a = aggregate(Expression_Rep1 ~ Gene, data = exp_a, FUN = max)
    # exp_b = exp2[exp2$Gene %in% enhancer$gene_name,]
    # exp_b = aggregate(Expression_Rep2 ~ Gene, data = exp_b, FUN = max)
    # dko1_exp = merge(exp_a, exp_b)

    DKO1_exp_list[[i]] <- exp
    corr <- cor(exp$Expression_Rep1, exp$Expression_Rep2, method = "spearman")
    DKO1_exp_corr_list[[i]] <- corr
  }
}

rm(exp1, exp2, exp_a, exp_b, exp, genanno, project_info, rse_gene, GROI)
names(DKO1_exp_list) <- c("cpgrich_MU", "cpgrich_UU", "cpgpoor_MU", "cpgpoor_UU") # "enhancer_MU", "enhancer_UU"
names(DKO1_exp_corr_list) <- c("cpgrich_MU", "cpgrich_UU", "cpgpoor_MU", "cpgpoor_UU") # "enhancer_MU", "enhancer_UU"

#-------------------------------------------------------------------------------
# HCT116 summarization function
#-------------------------------------------------------------------------------

HCT116_summarizeOverlaps = function(refGR, testGR, class, cutoff){
  
  #1 Find overlapping regions meeting a high cutoff
  hits <- findOverlaps(query = refGR, subject = testGR)
  overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
  hits <- hits[percentOverlap > cutoff]
  
  #2 Summarize overlaps to Ensembl id level
  df = data.frame(refGR$gene_id[queryHits(hits)], 
                  round(data.frame(mcols(testGR)[subjectHits(hits),]), 3), 
                  stringsAsFactors = F)
  
  if(ncol(df) > 2 & (class == "Enhancers" | class == "Promoters")){
    colnames(df)[1] = "Gene_id"
    df = t(sapply(split(df, df$Gene_id), function(x)x[which.max(x$HCT116_GSE126215_oATAC_Rep1.bw),]))
    df = as.data.frame(apply(df, 2, unlist))
    df[,2:ncol(df)] = apply(df[,2:ncol(df)], 2, as.numeric)
  }
  
  if(ncol(df) == 2){
    colnames(df) = c("Gene_id", "Methylation")
    df = aggregate(Methylation ~ Gene_id, data = df, FUN = median)
  }
  
  df = na.omit(df)
  return(df)
}

#-------------------------------------------------------------------------------
# HCT116 consolidated epigenomic matrix
#-------------------------------------------------------------------------------

HCT116_consolidateEpiMat = function(regRegion, summVals, regType, exp){
  
  hct116 = HCT116_summarizeOverlaps(refGR = regRegion, testGR = summVals, class = regType, cutoff = 0.999)
  hct116_meth = HCT116_summarizeOverlaps(refGR = regRegion, testGR = HCT116_meth, class = regType, cutoff = 0.999)
  
  hct116 = merge(hct116, hct116_meth)
  hct116 = merge(hct116, exp)
  rownames(hct116) = hct116$Gene_id # hct116$Gene for enhancer data
  hct116 = hct116[, -1]
  
  # Replicate1
  hct116_Rep1 = hct116[,c("Methylation", grep("Rep1", colnames(hct116), value = T))]
  hct116_Rep1 = hct116_Rep1[,c(10, 1,2,3:9)]
  colnames(hct116_Rep1) = c("Expression", "Methylation", "ATAC", "H2AZ", "H3K27ac", 
                            "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
  
  hct116_Rep1[hct116_Rep1 < 0] = 0
  
  hct116_Rep1 = log2(hct116_Rep1 + 1)
  hct116_Rep1 = na.omit(hct116_Rep1)
  hct116_Rep1 = hct116_Rep1[rowSums(hct116_Rep1) > 0,]
  
  # Replicate2
  hct116_Rep2 = hct116[,c("Methylation", grep("Rep2", colnames(hct116), value = T))]
  hct116_Rep2 = hct116_Rep2[,c(10, 1,2,3:9)]
  colnames(hct116_Rep2) = c("Expression", "Methylation", "ATAC", "H2AZ", "H3K27ac", 
                            "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
  
  hct116_Rep2[hct116_Rep2 < 0] = 0
  
  hct116_Rep2 = log2(hct116_Rep2 + 1)
  hct116_Rep2 = na.omit(hct116_Rep2)
  hct116_Rep2 = hct116_Rep2[is.finite(rowSums(hct116_Rep2)),]
  
  return(list(Replicate1 = hct116_Rep1, Replicate2 = hct116_Rep2))
}

#-------------------------------------------------------------------------------
# DKO1 summarization function
#-------------------------------------------------------------------------------

DKO1_summarizeOverlaps = function(refGR, testGR, class, cutoff){
  
  #1 Find overlapping regions meeting a high cutoff
  hits <- findOverlaps(query = refGR, subject = testGR)
  overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
  hits <- hits[percentOverlap > cutoff]
  
  #2 Summarize overlaps to Ensembl id level
  df = data.frame(refGR$gene_id[queryHits(hits)], 
                  round(data.frame(mcols(testGR)[subjectHits(hits),]), 3), 
                  stringsAsFactors = F)
  
  if(ncol(df) > 2 & (class == "Enhancers" | class == "Promoters")){
    colnames(df)[1] = "Gene_id"
    df = t(sapply(split(df, df$Gene_id), function(x)x[which.max(x$DKO1_GSE126215_oATAC_Rep1.bw),]))
    df = as.data.frame(apply(df, 2, unlist))
    df[,2:ncol(df)] = apply(df[,2:ncol(df)], 2, as.numeric)
  }
  
  if(ncol(df) == 2){
    colnames(df) = c("Gene_id", "Methylation")
    df = aggregate(Methylation ~ Gene_id, data = df, FUN = median)
  }
  
  df = na.omit(df)
  return(df)
}

#-------------------------------------------------------------------------------
# DKO1 consolidated epigenomic matrix
#-------------------------------------------------------------------------------

DKO1_consolidateEpiMat = function(regRegion, summVals, regType, exp){
  
  dko1 = DKO1_summarizeOverlaps(refGR = regRegion, testGR = summVals, class = regType, cutoff = 0.999)
  dko1_meth = DKO1_summarizeOverlaps(refGR = regRegion, testGR = DKO1_meth, class = regType, cutoff = 0.999)
  
  dko1 = merge(dko1, dko1_meth)
  dko1 = merge(dko1, exp)
  
  # For DKO1 cpgpoor, there are 4 data points having identical gene id
  if(length(dko1$Gene_id) != length(unique(dko1$Gene_id))){
    gene_id_count <- as.data.frame(table(dko1$Gene_id))
    gene_id_count <- gene_id_count[gene_id_count$Freq == 1,]
    dko1 <- dko1[dko1$Gene_id %in% gene_id_count$Var1,]
  }
  
  rownames(dko1) = dko1$Gene_id # dko1$Gene for enhancer data
  dko1 = dko1[, -1]
  
  # Replicate1
  dko1_Rep1 = dko1[,c("Methylation", grep("Rep1", colnames(dko1), value = T))]
  dko1_Rep1 = dko1_Rep1[,c(10, 1,2,3:9)]
  colnames(dko1_Rep1) = c("Expression", "Methylation", "ATAC", "H2AZ", "H3K27ac",
                          "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
  
  
  # For comparing Bigwig expression and Recount expression
  # dko1_Rep1 = dko1_Rep1[,c(11, 1,2,3:10)]
  # colnames(dko1_Rep1) = c("Expression", "Methylation", "ATAC", "BW_expression", "H2AZ", "H3K27ac", 
  #                         "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
  
  dko1_Rep1[dko1_Rep1 < 0] = 0
  
  dko1_Rep1 = log2(dko1_Rep1+1)
  dko1_Rep1 = na.omit(dko1_Rep1)
  dko1_Rep1 = dko1_Rep1[is.finite(rowSums(dko1_Rep1)),]
  
  # Replicate2
  dko1_Rep2 = dko1[,c("Methylation", grep("Rep2", colnames(dko1), value = T))]
  dko1_Rep2 = dko1_Rep2[,c(10, 1,2,3:9)]
  colnames(dko1_Rep2) = c("Expression", "Methylation", "ATAC", "H2AZ", "H3K27ac",
                          "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
  
  
  # For comparing Bigwig expression and Recount expression
  # dko1_Rep1 = dko1_Rep1[,c(11, 1,2,3:10)]
  # colnames(dko1_Rep1) = c("Expression", "Methylation", "ATAC", "BW_expression", "H2AZ", "H3K27ac", 
  #                         "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
  
  dko1_Rep2[dko1_Rep2 < 0] = 0
  
  dko1_Rep2 = log2(dko1_Rep2+1)
  dko1_Rep2 = na.omit(dko1_Rep2)
  dko1_Rep2 = dko1_Rep2[is.finite(rowSums(dko1_Rep2)),]
  
  return(list(Replicate1 = dko1_Rep1, Replicate2 = dko1_Rep2))
}

#-------------------------------------------------------------------------------
# Run
#-------------------------------------------------------------------------------

hct116_cpgrich_MU = HCT116_consolidateEpiMat(regRegion = cpgrich_MU, summVals = HCT116_cpgrich_MU, regType = "Promoters", exp = HCT116_exp_list[[1]])
hct116_cpgrich_UU = HCT116_consolidateEpiMat(regRegion = cpgrich_UU, summVals = HCT116_cpgrich_UU, regType = "Promoters", exp = HCT116_exp_list[[2]])
#hct116_enhancer_MU = HCT116_consolidateEpiMat(regRegion = enhancer_MU, summVals = HCT116_enhancer_MU, regType = "Enhancers", exp = HCT116_exp_list[[3]])
#hct116_enhancer_UU = HCT116_consolidateEpiMat(regRegion = enhancer_UU, summVals = HCT116_enhancer_UU, regType = "Enhancers", exp = HCT116_exp_list[[4]])
hct116_cpgpoor_MU = HCT116_consolidateEpiMat(regRegion = cpgpoor_MU, summVals = HCT116_cpgpoor_MU, regType = "Promoters", exp = HCT116_exp_list[[2]])
#hct116_cpgpoor_UU = HCT116_consolidateEpiMat(regRegion = cpgpoor_UU, summVals = HCT116_cpgpoor_UU, regType = "Promoters", exp = HCT116_exp_list[[4]])
dko1_cpgrich_MU = DKO1_consolidateEpiMat(regRegion = cpgrich_MU, summVals = DKO1_cpgrich_MU, regType = "Promoters", exp = DKO1_exp_list[[1]])
dko1_cpgrich_UU = DKO1_consolidateEpiMat(regRegion = cpgrich_UU, summVals = DKO1_cpgrich_UU, regType = "Promoters", exp = DKO1_exp_list[[2]])
#dko1_enhancer_MU = DKO1_consolidateEpiMat(regRegion = enhancer_MU, summVals = DKO1_enhancer_MU, regType = "Enhancers", exp = DKO1_exp_list[[3]])
#dko1_enhancer_UU = DKO1_consolidateEpiMat(regRegion = enhancer_UU, summVals = DKO1_enhancer_UU, regType = "Enhancers", exp = DKO1_exp_list[[4]])
dko1_cpgpoor_MU = DKO1_consolidateEpiMat(regRegion = cpgpoor_MU, summVals = DKO1_cpgpoor_MU, regType = "Promoters", exp = DKO1_exp_list[[2]])
dko1_cpgpoor_UU = DKO1_consolidateEpiMat(regRegion = cpgpoor_UU, summVals = DKO1_cpgpoor_UU, regType = "Promoters", exp = DKO1_exp_list[[4]])

#-------------------------------------------------------------------------------
# Exclude regions having improper methylation values
#-------------------------------------------------------------------------------

hct116_cpgrich_MU <- list(Replicate1 = hct116_cpgrich_MU$Replicate1[which(hct116_cpgrich_MU$Replicate1$Methylation >= 5),],
                          Replicate2 = hct116_cpgrich_MU$Replicate2[which(hct116_cpgrich_MU$Replicate2$Methylation >= 5),])
hct116_cpgrich_UU <- list(Replicate1 = hct116_cpgrich_UU$Replicate1[which(hct116_cpgrich_UU$Replicate1$Methylation < 5),],
                          Replicate2 = hct116_cpgrich_UU$Replicate2[which(hct116_cpgrich_UU$Replicate2$Methylation < 5),])
hct116_cpgpoor_MU <- list(Replicate1 = hct116_cpgpoor_MU$Replicate1[which(hct116_cpgpoor_MU$Replicate1$Methylation >= 5),],
                         Replicate2 = hct116_cpgpoor_MU$Replicate2[which(hct116_cpgpoor_MU$Replicate2$Methylation >= 5),])
hct116_cpgpoor_UU <- list(Replicate1 = hct116_cpgpoor_UU$Replicate1[which(hct116_cpgpoor_UU$Replicate1$Methylation < 5),],
                         Replicate2 = hct116_cpgpoor_UU$Replicate2[which(hct116_cpgpoor_UU$Replicate2$Methylation < 5),])
dko1_cpgrich_MU <- list(Replicate1 = dko1_cpgrich_MU$Replicate1[which(dko1_cpgrich_MU$Replicate1$Methylation < 5),],
                        Replicate2 = dko1_cpgrich_MU$Replicate2[which(dko1_cpgrich_MU$Replicate2$Methylation < 5),])
dko1_cpgrich_UU <- list(Replicate1 = dko1_cpgrich_UU$Replicate1[which(dko1_cpgrich_UU$Replicate1$Methylation < 5),],
                        Replicate2 = dko1_cpgrich_UU$Replicate2[which(dko1_cpgrich_UU$Replicate2$Methylation < 5),])
dko1_cpgpoor_MU <- list(Replicate1 = dko1_cpgpoor_MU$Replicate1[which(dko1_cpgpoor_MU$Replicate1$Methylation < 5),],
                       Replicate2 = dko1_cpgpoor_MU$Replicate2[which(dko1_cpgpoor_MU$Replicate2$Methylation < 5),])
dko1_cpgpoor_UU <- list(Replicate1 = dko1_cpgpoor_UU$Replicate1[which(dko1_cpgpoor_UU$Replicate1$Methylation < 5),],
                       Replicate2 = dko1_cpgpoor_UU$Replicate2[which(dko1_cpgpoor_UU$Replicate2$Methylation < 5),])

#-------------------------------------------------------------------------------
# Sanity check if Rep1 and Rep2 share the same regions
#-------------------------------------------------------------------------------

identical(rownames(hct116_cpgrich_MU$Replicate1), rownames(hct116_cpgrich_MU$Replicate2))
identical(rownames(hct116_cpgrich_UU$Replicate1), rownames(hct116_cpgrich_UU$Replicate2))
identical(rownames(hct116_cpgpoor_MU$Replicate1), rownames(hct116_cpgpoor_MU$Replicate2))
identical(rownames(hct116_cpgpoor_UU$Replicate1), rownames(hct116_cpgpoor_UU$Replicate2))
identical(rownames(dko1_cpgrich_MU$Replicate1), rownames(dko1_cpgrich_MU$Replicate2))
identical(rownames(dko1_cpgrich_UU$Replicate1), rownames(dko1_cpgrich_UU$Replicate2))
identical(rownames(dko1_cpgpoor_MU$Replicate1), rownames(dko1_cpgpoor_MU$Replicate2))
identical(rownames(dko1_cpgpoor_UU$Replicate1), rownames(dko1_cpgpoor_UU$Replicate2))

#-------------------------------------------------------------------------------
# Pick out Common regions between HCT116 and DKO1
#-------------------------------------------------------------------------------

cmn_hct116 <- rownames(hct116_cpgrich_MU$Replicate1) %in% rownames(dko1_cpgrich_MU$Replicate1)
cmn_dko1 <- rownames(dko1_cpgrich_MU$Replicate1) %in% rownames(hct116_cpgrich_MU$Replicate1)
hct116_cpgrich_MU <- list(Replicate1 = hct116_cpgrich_MU$Replicate1[cmn_hct116 == T,], Replicate2 = hct116_cpgrich_MU$Replicate2[cmn_hct116 == T,])
dko1_cpgrich_MU <- list(Replicate1 = dko1_cpgrich_MU$Replicate1[cmn_dko1 == T,], Replicate2 = dko1_cpgrich_MU$Replicate2[cmn_dko1 == T,])

cmn_hct116 <- rownames(hct116_cpgrich_UU$Replicate1) %in% rownames(dko1_cpgrich_UU$Replicate1)
cmn_dko1 <- rownames(dko1_cpgrich_UU$Replicate1) %in% rownames(hct116_cpgrich_UU$Replicate1)
hct116_cpgrich_UU <- list(Replicate1 = hct116_cpgrich_UU$Replicate1[cmn_hct116 == T,], Replicate2 = hct116_cpgrich_UU$Replicate2[cmn_hct116 == T,])
dko1_cpgrich_UU <- list(Replicate1 = dko1_cpgrich_UU$Replicate1[cmn_dko1 == T,], Replicate2 = dko1_cpgrich_UU$Replicate2[cmn_dko1 == T,])

cmn_hct116 <- rownames(hct116_cpgpoor_MU$Replicate1) %in% rownames(dko1_cpgpoor_MU$Replicate1)
cmn_dko1 <- rownames(dko1_cpgpoor_MU$Replicate1) %in% rownames(hct116_cpgpoor_MU$Replicate1)
hct116_cpgpoor_MU <- list(Replicate1 = hct116_cpgpoor_MU$Replicate1[cmn_hct116 == T,], Replicate2 = hct116_cpgpoor_MU$Replicate2[cmn_hct116 == T,])
dko1_cpgpoor_MU <- list(Replicate1 = dko1_cpgpoor_MU$Replicate1[cmn_dko1 == T,], Replicate2 = dko1_cpgpoor_MU$Replicate2[cmn_dko1 == T,])

cmn_hct116 <- rownames(hct116_cpgpoor_UU$Replicate1) %in% rownames(dko1_cpgpoor_UU$Replicate1)
cmn_dko1 <- rownames(dko1_cpgpoor_UU$Replicate1) %in% rownames(hct116_cpgpoor_UU$Replicate1)
hct116_cpgpoor_UU <- list(Replicate1 = hct116_cpgpoor_UU$Replicate1[cmn_hct116 == T,], Replicate2 = hct116_cpgpoor_UU$Replicate2[cmn_hct116 == T,])
dko1_cpgpoor_UU <- list(Replicate1 = dko1_cpgpoor_UU$Replicate1[cmn_dko1 == T,], Replicate2 = dko1_cpgpoor_UU$Replicate2[cmn_dko1 == T,])

#-------------------------------------------------------------------------------
# Sanity check if HCT116 and DKO1 share the same regions
#-------------------------------------------------------------------------------

identical(dim(hct116_cpgrich_MU$Replicate1), dim(dko1_cpgrich_MU$Replicate1))
identical(dim(hct116_cpgrich_UU$Replicate1), dim(dko1_cpgrich_UU$Replicate1))
identical(dim(hct116_cpgpoor_MU$Replicate1), dim(dko1_cpgpoor_MU$Replicate1))
identical(dim(hct116_cpgpoor_UU$Replicate1), dim(dko1_cpgpoor_UU$Replicate1))

#-------------------------------------------------------------------------------
# Save data
#-------------------------------------------------------------------------------

if(!dir.exists(paste0(path2, "summarizedMatrix/All_Consolidated_HCT116/"))){
  dir.create(paste0(path2, "summarizedMatrix/All_Consolidated_HCT116/"), recursive = T)
}
if(!dir.exists(paste0(path2, "summarizedMatrix/All_Consolidated_DKO1/"))){
  dir.create(paste0(path2, "summarizedMatrix/All_Consolidated_DKO1/"), recursive = T)
}

saveRDS(hct116_cpgrich_MU, paste0(path2, "summarizedMatrix/All_Consolidated_HCT116/hct116_cpgrich_MU_allEpigenomes.RDS"))
saveRDS(hct116_cpgrich_UU, paste0(path2, "summarizedMatrix/All_Consolidated_HCT116/hct116_cpgrich_UU_allEpigenomes.RDS"))
#saveRDS(hct116_enhancer_MU, paste0(path2, "summarizedMatrix/All_Consolidated_HCT116/hct116_enhancer_MU_allEpigenomes.RDS"))
#saveRDS(hct116_enhancer_UU, paste0(path2, "summarizedMatrix/All_Consolidated_HCT116/hct116_enhancer_UU_allEpigenomes.RDS"))
saveRDS(hct116_cpgpoor_MU, paste0(path2, "summarizedMatrix/All_Consolidated_HCT116/hct116_cpgpoor_MU_allEpigenomes.RDS"))
saveRDS(hct116_cpgpoor_UU, paste0(path2, "summarizedMatrix/All_Consolidated_HCT116/hct116_cpgpoor_UU_allEpigenomes.RDS"))
saveRDS(dko1_cpgrich_MU, paste0(path2, "summarizedMatrix/All_Consolidated_DKO1/dko1_cpgrich_MU_allEpigenomes.RDS"))
saveRDS(dko1_cpgrich_UU, paste0(path2, "summarizedMatrix/All_Consolidated_DKO1/dko1_cpgrich_UU_allEpigenomes.RDS"))
#saveRDS(dko1_enhancer_MU, paste0(path2, "summarizedMatrix/All_Consolidated_DKO1/dko1_enhancer_MU_allEpigenomes.RDS"))
#saveRDS(dko1_enhancer_UU, paste0(path2, "summarizedMatrix/All_Consolidated_DKO1/dko1_enhancer_UU_allEpigenomes.RDS"))
saveRDS(dko1_cpgpoor_MU, paste0(path2, "summarizedMatrix/All_Consolidated_DKO1/dko1_cpgpoor_MU_allEpigenomes.RDS"))
saveRDS(dko1_cpgpoor_UU, paste0(path2, "summarizedMatrix/All_Consolidated_DKO1/dko1_cpgpoor_UU_allEpigenomes.RDS"))

################################################################################
# Extract regions in specific consolidated matrix from BED file
################################################################################

# Sanity check
# identical(rownames(hct116_cpgrich_MU$Replicate1), rownames(hct116_cpgrich_MU$Replicate2))
# identical(rownames(hct116_cpgrich_MU$Replicate1), rownames(dko1_cpgrich_MU$Replicate1))
# identical(rownames(hct116_cpgrich_UU$Replicate1), rownames(hct116_cpgrich_UU$Replicate2))
# identical(rownames(hct116_cpgrich_UU$Replicate1), rownames(dko1_cpgrich_UU$Replicate1))
# identical(rownames(hct116_cpgrich$Replicate1), rownames(hct116_cpgrich$Replicate2))
# identical(rownames(hct116_cpgpoor$Replicate1), rownames(hct116_cpgpoor$Replicate1))
# identical(rownames(dko1_cpgrich$Replicate1), rownames(dko1_cpgrich$Replicate2))
# identical(rownames(dko1_cpgpoor$Replicate1), rownames(dko1_cpgpoor$Replicate1))

# HCT116 CGI promoter MU
MU <- vector()
for(i in 1:nrow(hct116_cpgrich_MU$Replicate1)){
  n <- which(mcols(cpgrich_MU)$gene_id == rownames(hct116_cpgrich_MU$Replicate1)[i])
  MU[i] <- n
}
MU <- sort(MU)
HCT116_cpgrich_MU <- cpgrich_MU[MU]

# HCT116 CGI promoter UU
UU <- vector()
for(i in 1:nrow(hct116_cpgrich_UU$Replicate1)){
  n <- which(mcols(cpgrich_UU)$gene_id == rownames(hct116_cpgrich_UU$Replicate1)[i])
  UU[i] <- n
}
UU <- sort(UU)
HCT116_cpgrich_UU <- cpgrich_UU[UU]

# DKO1 CGI promoter MU
MU <- vector()
for(i in 1:nrow(dko1_cpgrich_MU$Replicate1)){
  n <- which(mcols(cpgrich_MU)$gene_id == rownames(dko1_cpgrich_MU$Replicate1)[i])
  MU[i] <- n
}
MU <- sort(MU)
DKO1_cpgrich_MU <- cpgrich_MU[MU]

# DKO1 CGI promoter UU
UU <- vector()
for(i in 1:nrow(dko1_cpgrich_UU$Replicate1)){
  n <- which(mcols(cpgrich_UU)$gene_id == rownames(dko1_cpgrich_UU$Replicate1)[i])
  UU[i] <- n
}
UU <- sort(UU)
DKO1_cpgrich_UU <- cpgrich_UU[UU]

# HCT116 non-CGI promoter MU
MU <- vector()
for(i in 1:nrow(hct116_cpgpoor_MU$Replicate1)){
  n <- which(mcols(cpgpoor_MU)$gene_id == rownames(hct116_cpgpoor_MU$Replicate1)[i])
  MU[i] <- n
}
MU <- sort(MU)
HCT116_cpgpoor_MU <- cpgpoor_MU[MU]

# HCT116 non-CGI promoter UU
UU <- vector()
for(i in 1:nrow(hct116_cpgpoor_UU$Replicate1)){
  n <- which(mcols(cpgpoor_UU)$gene_id == rownames(hct116_cpgpoor_UU$Replicate1)[i])
  UU[i] <- n
}
UU <- sort(UU)
HCT116_cpgpoor_UU <- cpgpoor_UU[UU]

# DKO1 non-CGI promoter MU
MU <- vector()
for(i in 1:nrow(dko1_cpgpoor_MU$Replicate1)){
  n <- which(mcols(cpgpoor_MU)$gene_id == rownames(dko1_cpgpoor_MU$Replicate1)[i])
  MU[i] <- n
}
MU <- sort(MU)
DKO1_cpgpoor_MU <- cpgpoor_MU[MU]

# DKO1 non-CGI promoter UU
UU <- vector()
for(i in 1:nrow(dko1_cpgpoor_UU$Replicate1)){
  n <- which(mcols(cpgpoor_UU)$gene_id == rownames(dko1_cpgpoor_UU$Replicate1)[i])
  UU[i] <- n
}
UU <- sort(UU)
DKO1_cpgpoor_UU <- cpgpoor_UU[UU]

################################################################################
# Save data
# Constitute of genomic regions of HCT116 equals to DKO1 if they share the same regulatory regions
# (depends on "Pick out Common regions between HCT116 and DKO1")
################################################################################

path4 <- "/path/to/dir/" # Save results for ChIPseeker visualization

if(!dir.exists(path4)){
  dir.create(path4, recursive = T)
}

export.bed(object = HCT116_cpgrich_MU, con = paste0(path4, "HCT116_cpgrich_MU.bed"))
saveRDS(HCT116_cpgrich_MU, paste0(path4, "HCT116_cpgrich_MU.RDS"))

export.bed(object = HCT116_cpgrich_UU, con = paste0(path4, "HCT116_cpgrich_UU.bed"))
saveRDS(HCT116_cpgrich_UU, paste0(path4, "HCT116_cpgrich_UU.RDS"))

export.bed(object = DKO1_cpgrich_MU, con = paste0(path4, "DKO1_cpgrich_MU.bed"))
saveRDS(DKO1_cpgrich_MU, paste0(path4, "DKO1_cpgrich_MU.RDS"))

export.bed(object = DKO1_cpgrich_UU, con = paste0(path4, "DKO1_cpgrich_UU.bed"))
saveRDS(DKO1_cpgrich_UU, paste0(path4, "DKO1_cpgrich_UU.RDS"))

export.bed(object = HCT116_cpgpoor_MU, con = paste0(path4, "HCT116_cpgpoor_MU.bed"))
saveRDS(HCT116_cpgpoor_MU, paste0(path4, "HCT116_cpgpoor_MU.RDS"))

export.bed(object = HCT116_cpgpoor_UU, con = paste0(path4, "HCT116_cpgpoor_UU.bed"))
saveRDS(HCT116_cpgpoor_UU, paste0(path4, "HCT116_cpgpoor_UU.RDS"))

export.bed(object = DKO1_cpgpoor_MU, con = paste0(path4, "DKO1_cpgpoor_MU.bed"))
saveRDS(DKO1_cpgpoor_MU, paste0(path4, "DKO1_cpgpoor_MU.RDS"))

export.bed(object = DKO1_cpgpoor_UU, con = paste0(path4, "DKO1_cpgpoor_UU.bed"))
saveRDS(DKO1_cpgpoor_UU, paste0(path4, "DKO1_cpgpoor_UU.RDS"))
