path <- "/path/to/summarizedMatrix/"
data_path <- "/path/to/dir/"

# Create output directory
if(!dir.exists(paste0(path, "plotMethAcceHeatmap"))){
  dir.create(paste0(path, "plotMethAcceHeatmap"), recursive = T)}

################################################################################
# Generate computed matrix of promoters over DNA methylation and accessibility
# for plotHeatmap
################################################################################

# CpG-rich MU
system(paste("computeMatrix", "scale-regions",
             "--scoreFileName",
             paste0(data_path, "NOMe_seq/HCT116_NOMe_HCG.bw"),
             paste0(data_path, "NOMe_seq/DKO1_NOMe_HCG.bw"),
             paste0(data_path, "NOMe_seq/HCT116_NOMe_GCH.bw"),
             paste0(data_path, "NOMe_seq/DKO1_NOMe_GCH.bw"),
             "--regionsFileName", paste0(data_path, "BED/ChIPseeker/cpgrich_MU.bed"),
             "--regionBodyLength", 3000, "--binSize", 10,
             "--outFileName", paste0(path, "plotMethAcceHeatmap/cpgrich_MU_matrix.npz"),
             "--outFileNameMatrix", paste0(path, "plotMethAcceHeatmap/cpgrich_MU_matrix.tab")))
print("CpG-rich MU matrix done")

# CpG-rich UU
system(paste("computeMatrix", "scale-regions",
             "--scoreFileName",
             paste0(data_path, "NOMe_seq/HCT116_NOMe_HCG.bw"),
             paste0(data_path, "NOMe_seq/DKO1_NOMe_HCG.bw"),
             paste0(data_path, "NOMe_seq/HCT116_NOMe_GCH.bw"),
             paste0(data_path, "NOMe_seq/DKO1_NOMe_GCH.bw"),
             "--regionsFileName", paste0(data_path, "BED/ChIPseeker/cpgrich_UU.bed"),
             "--regionBodyLength", 3000, "--binSize", 10,
             "--outFileName", paste0(path, "plotMethAcceHeatmap/cpgrich_UU_matrix.npz"),
             "--outFileNameMatrix", paste0(path, "plotMethAcceHeatmap/cpgrich_UU_matrix.tab")))
print("CpG-rich UU matrix done")

# CpG-poor MU
system(paste("computeMatrix", "scale-regions",
             "--scoreFileName",
             paste0(data_path, "NOMe_seq/HCT116_NOMe_HCG.bw"),
             paste0(data_path, "NOMe_seq/DKO1_NOMe_HCG.bw"),
             paste0(data_path, "NOMe_seq/HCT116_NOMe_GCH.bw"),
             paste0(data_path, "NOMe_seq/DKO1_NOMe_GCH.bw"),
             "--regionsFileName", paste0(data_path, "BED/ChIPseeker/cpgpoor_MU.bed"),
             "--regionBodyLength", 3000, "--binSize", 10,
             "--outFileName", paste0(path, "plotMethAcceHeatmap/cpgpoor_MU_matrix.npz"),
             "--outFileNameMatrix", paste0(path, "plotMethAcceHeatmap/cpgpoor_MU_matrix.tab")))
print("CpG-poor MU matrix done")

################################################################################
# Plot enriched heatmap for DNA methylation and accessibility
################################################################################

# CpG-rich MU
system(paste("plotHeatmap",
             "--matrixFile", paste0(path, "plotMethAcceHeatmap/cpgrich_MU_matrix.npz"),
             "--outFileName", paste0(path, "plotMethAcceHeatmap/cpgrich_MU_heatmap.png"),
             "--interpolationMethod", "nearest",
             "--colorList", "black,yellow",
             "--heatmapHeight", 16, "--heatmapWidth", 8,
             "--boxAroundHeatmaps", "no",
             "--xAxisLabel", "'Distance to TSS (bp)'",
             "--startLabel", "'-1500'", "--endLabel", "'1500'",
             "--regionsLabel", "'Genes'",
             "--legendLocation", "none",
             "--samplesLabel", "'HCT116-DNA Methylation'", "'DKO1-DNA Methylation'",
             "'HCT116-Accessibility'", "'DKO1-Accessibility'",
             "--plotTitle", "'CGI promoter MU'"))
print("CpG-rich MU heatmap done")

# CpG-rich UU
system(paste("plotHeatmap",
             "--matrixFile", paste0(path, "plotMethAcceHeatmap/cpgrich_UU_matrix.npz"),
             "--outFileName", paste0(path, "plotMethAcceHeatmap/cpgrich_UU_heatmap.png"),
             "--interpolationMethod", "nearest",
             "--colorList", "black,yellow",
             "--heatmapHeight", 16, "--heatmapWidth", 8,
             "--boxAroundHeatmaps", "no",
             "--xAxisLabel", "'Distance to TSS (bp)'",
             "--startLabel", "'-1500'", "--endLabel", "'1500'",
             "--regionsLabel", "'Genes'",
             "--legendLocation", "none",
             "--samplesLabel", "'HCT116-DNA Methylation'", "'DKO1-DNA Methylation'",
             "'HCT116-Accessibility'", "'DKO1-Accessibility'",
             "--plotTitle", "'CGI promoter UU'"))
print("CpG-rich UU heatmap done")

# CpG-poor MU
system(paste("plotHeatmap",
             "--matrixFile", paste0(path, "plotMethAcceHeatmap/cpgpoor_MU_matrix.npz"),
             "--outFileName", paste0(path, "plotMethAcceHeatmap/cpgpoor_MU_heatmap.png"),
             "--interpolationMethod", "nearest",
             "--colorList", "black,yellow",
             "--heatmapHeight", 16, "--heatmapWidth", 8,
             "--boxAroundHeatmaps", "no",
             "--xAxisLabel", "'Distance to TSS (bp)'",
             "--startLabel", "'-1500'", "--endLabel", "'1500'",
             "--regionsLabel", "'Genes'",
             "--legendLocation", "none",
             "--samplesLabel", "'HCT116-DNA Methylation'", "'DKO1-DNA Methylation'",
             "'HCT116-Accessibility'", "'DKO1-Accessibility'",
             "--plotTitle", "'Non-CGI promoter MU'"))
print("CpG-poor MU heatmap done")

