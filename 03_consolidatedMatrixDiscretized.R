library(Ckmeans.1d.dp)
library(bnlearn)

path = "/path/to/dir/"

# Define regions of interest
cellType <- c("HCT116", "DKO1")
regType <- c("hct116_cpgrich_MU", "hct116_cpgrich_UU", "hct116_cpgpoor_MU", "hct116_cpgpoor_UU",
             "dko1_cpgrich_MU", "dko1_cpgrich_UU", "dko1_cpgpoor_MU", "dko1_cpgpoor_UU")

# Assign varied numbers of clusters for features
# dis_k <- matrix(nrow = 1, ncol = 10)
# rownames(dis_k) <- "k"
# colnames(dis_k) <- c("Expression", "Methylation", "ATAC", "H2AZ", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
# dis_k[1, 1:10] <- c(6, 3, 2, 6, 6, 6, 6, 6, 6, 6)

################################################################################
# Ckmeans.1d.dp
################################################################################

discretizeData = function(cellType, regType){
  # Access to consolidated matrices
  data_path = paste0(path, "summarizedMatrix/All_Consolidated_", cellType, "/", regType, "_allEpigenomes.RDS")
  dat = readRDS(data_path)
  
  dis = lapply(dat, function(val) {
    tmpdat = apply(val, 2, function(m) Ckmeans.1d.dp(x = m, k = 3)$cluster)
    rownames(tmpdat) = rownames(val)
    return(tmpdat)
  })

  # r1 <- dat$Replicate1
  # r2 <- dat$Replicate2
  # disDat_1 <- matrix(nrow = nrow(r1), ncol = ncol(r1))
  # rownames(disDat_1) <- rownames(r1)
  # colnames(disDat_1) <- colnames(r1)
  # disDat_2 <- matrix(nrow = nrow(r2), ncol = ncol(r2))
  # rownames(disDat_2) <- rownames(r2)
  # colnames(disDat_2) <- colnames(r2)
  
  # for(i in 1:length(dis_k)){
  #   tmp_1 <- Ckmeans.1d.dp(r1[,i], k = dis_k[i])$cluster
  #   disDat_1[,i] <- tmp_1
    
  #   tmp_2 <- Ckmeans.1d.dp(r2[,i], k = dis_k[i])$cluster
  #   disDat_2[,i] <- tmp_2
  # }
  # dis <- list(Replicate1 = disDat_1, Replicate2 = disDat_2)

  return(dis)
}

# Run
hct116_cpgrich_MU = discretizeData(cellType = cellType[1], regType = regType[1])
hct116_cpgrich_UU = discretizeData(cellType = cellType[1], regType = regType[2])
# hct116_enhancer_MU = discretizeData(cellType = cellType[1], regType = regType[3])
# hct116_enhancer_UU = discretizeData(cellType = cellType[1], regType = regType[4])
hct116_cpgpoor_MU = discretizeData(cellType = cellType[1], regType = regType[3])
hct116_cpgpoor_UU = discretizeData(cellType = cellType[1], regType = regType[4])
dko1_cpgrich_MU = discretizeData(cellType = cellType[2], regType = regType[5])
dko1_cpgrich_UU = discretizeData(cellType = cellType[2], regType = regType[6])
# dko1_enhancer_MU = discretizeData(cellType = cellType[2], regType = regType[11])
# dko1_enhancer_UU = discretizeData(cellType = cellType[2], regType = regType[12])
dko1_cpgpoor_MU = discretizeData(cellType = cellType[2], regType = regType[7])
dko1_cpgpoor_UU = discretizeData(cellType = cellType[2], regType = regType[8])

# Create the output dirs
if (!dir.exists(paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/"))) {
  dir.create(paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/"), recursive = T)
  }

# Save data
saveRDS(hct116_cpgrich_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/hct116_cpgrich_MU_discretized.RDS"))
saveRDS(hct116_cpgrich_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/hct116_cpgrich_UU_discretized.RDS"))
# saveRDS(hct116_enhancer_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/hct116_enhancer_MU_discretized.RDS"))
# saveRDS(hct116_enhancer_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/hct116_enhancer_UU_discretized.RDS"))
saveRDS(hct116_cpgpoor_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/hct116_cpgpoor_MU_discretized.RDS"))
saveRDS(hct116_cpgpoor_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/hct116_cpgpoor_UU_discretized.RDS"))
saveRDS(dko1_cpgrich_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/dko1_cpgrich_MU_discretized.RDS"))
saveRDS(dko1_cpgrich_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/dko1_cpgrich_UU_discretized.RDS"))
# saveRDS(dko1_enhancer_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/dko1_enhancer_MU_discretized.RDS"))
# saveRDS(dko1_enhancer_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/dko1_enhancer_UU_discretized.RDS"))
saveRDS(dko1_cpgpoor_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/dko1_cpgpoor_MU_discretized.RDS"))
saveRDS(dko1_cpgpoor_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/dko1_cpgpoor_UU_discretized.RDS"))

################################################################################
# Hartemink's method
################################################################################

discretizeDataHartemink <- function(cellType, regType){
  data_path <- paste0(path, "summarizedMatrix/All_Consolidated_", cellType, "/", regType, "_allEpigenomes.RDS")
  dat <- readRDS(data_path)
  
  dis <- lapply(dat, function(val){
    tmpdat =  discretize(val, method = "hartemink", breaks = 3, ibreaks = 60, idisc = "quantile")
    return(tmpdat)
    })
  return(dis)
}

# Run
hct116_cpgrich_MU = discretizeDataHartemink(cellType = cellType[1], regType = regType[1])
hct116_cpgrich_UU = discretizeDataHartemink(cellType = cellType[1], regType = regType[2])
hct116_cpgpoor_MU = discretizeDataHartemink(cellType = cellType[1], regType = regType[3])
hct116_cpgpoor_UU = discretizeDataHartemink(cellType = cellType[1], regType = regType[4])
dko1_cpgrich_MU = discretizeDataHartemink(cellType = cellType[2], regType = regType[5])
dko1_cpgrich_UU = discretizeDataHartemink(cellType = cellType[2], regType = regType[6])
dko1_cpgpoor_MU = discretizeDataHartemink(cellType = cellType[2], regType = regType[7])
dko1_cpgpoor_UU = discretizeDataHartemink(cellType = cellType[2], regType = regType[8])

# Create the output dirs
if (!dir.exists(paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/"))) {
  dir.create(paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/"), recursive = T)
}

# Save data
saveRDS(hct116_cpgrich_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/hct116_cpgrich_MU_discretized.RDS"))
saveRDS(hct116_cpgrich_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/hct116_cpgrich_UU_discretized.RDS"))
saveRDS(hct116_cpgpoor_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/hct116_cpgpoor_MU_discretized.RDS"))
saveRDS(hct116_cpgpoor_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/hct116_cpgpoor_UU_discretized.RDS"))
saveRDS(dko1_cpgrich_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/dko1_cpgrich_MU_discretized.RDS"))
saveRDS(dko1_cpgrich_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/dko1_cpgrich_UU_discretized.RDS"))
saveRDS(dko1_cpgpoor_MU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/dko1_cpgpoor_MU_discretized.RDS"))
saveRDS(dko1_cpgpoor_UU, paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized_Hartemink/dko1_cpgpoor_UU_discretized.RDS"))
