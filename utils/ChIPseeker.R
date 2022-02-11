library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

path <- "/path/to/dir/" #Access to BED files
files <- c("HCT116_cpgrich_MU.bed", "DKO1_cpgrich_MU.bed", "HCT116_cpgpoor_MU.bed", "DKO1_cpgpoor_MU.bed")

# Create output directory
if(!dir.exists("/path/to/ChIPseeker")){
  dir.create("/path/to/ChIPseeker")
}

################################################################################
# Annotation pie plot
################################################################################

for(i in 1:length(files)){
  peakAnno <- annotatePeak(paste0(path, files[i]), tssRegion = c(-1500, 1500), TxDb = txdb)
  plotAnnoPie(peakAnno)
}

################################################################################
# Profile of peaks corresponding to TSS regions
################################################################################

num_row_list <- list()
for(i in 1:length(files)){
  peakHeat <- peakHeatmap(paste0(path, files[i]), TxDb = txdb, upstream = 3000, downstream = 3000, color = "red", xlab = "Distance to TSS (bp)")
  num_row_list[[i]] <- nrow(peakHeat)
}
