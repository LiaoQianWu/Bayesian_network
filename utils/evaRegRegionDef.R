path <- "/path/to/summarizedMatrix/"

hct116_cpgrich_MU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgrich_MU_allEpigenomes.RDS"))
hct116_cpgrich_UU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgrich_UU_allEpigenomes.RDS"))
hct116_cpgpoor_MU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgpoor_MU_allEpigenomes.RDS"))
hct116_cpgpoor_UU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgpoor_UU_allEpigenomes.RDS"))
dko1_cpgrich_MU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgrich_MU_allEpigenomes.RDS"))
dko1_cpgrich_UU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgrich_UU_allEpigenomes.RDS"))
dko1_cpgpoor_MU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgpoor_MU_allEpigenomes.RDS"))
dko1_cpgpoor_UU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgpoor_UU_allEpigenomes.RDS"))

regType <- c(hct116_cpgrich_MU, hct116_cpgrich_UU,
             hct116_cpgpoor_MU, hct116_cpgpoor_UU,
             dko1_cpgrich_MU, dko1_cpgrich_UU,
             dko1_cpgpoor_MU, dko1_cpgpoor_UU)
# regType <- regType[seq(1, length(regType), 2)]

# Dim
for(i in regType){
  print(dim(i))
}

# Unmethylation rate
for(i in regType){
  unmeth <- length(which(i$Methylation < 5))
  meth_all <- length(i$Methylation)
  unmeth_ratio <- unmeth / meth_all
  print(unmeth_ratio)
}

# Overlapping sites
i <- 1
j <- 2
for(p in 1:4){
  num <- length(which(rownames(regType[[i]]) %in% rownames(regType[[j]]) == T))
  print(num)
  i <- i + 2
  j <- j + 2
}

