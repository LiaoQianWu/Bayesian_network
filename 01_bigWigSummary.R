# Initialize the arguments for multiple runs
args = commandArgs(trailingOnly = T)
i = as.numeric(args[1])
j = as.numeric(args[2])


path = "/path/to/dir/" #Access to BED files
path2 = "/path/to/dir/" #Save summarized results

# List data path for running deepTools commands
bw_imputed_path = list.dirs("/path/to/dir", recursive = F) #Access to e.g., ChIP-seq, RNA-seq, or WGBS data
bw_imputed_path = bw_imputed_path[basename(bw_imputed_path) != "data"]

bw_files = lapply(bw_imputed_path, function(samp){
  bw = list.files(samp, full.names = T)
  bw = grep(".bw", bw, value = T)
})
names(bw_files) = basename(bw_imputed_path)

reg_type = c(
  cpgrich_MU = paste0(path, "cpgrich_MU.bed"),
  cpgrich_UU = paste0(path, "cpgrich_UU.bed"),
  enhancer_MU = paste0(path, "enhancer_MU.bed"),
  enhancer_UU = paste0(path, "enhancer_UU.bed"),
  cpgpoor_MU = paste0(path, "cpgpoor_MU.bed"),
  cpgpoor_UU = paste0(path, "cpgpoor_UU.bed")
)

# Create output directory structure
if(!dir.exists(paste0(path2, "multiBigwigSummary/"))){
  dir.create(paste0(path2, "multiBigwigSummary/"), recursive = T)
  
  sapply(names(bw_files), function(x){
    dir.create(paste0(path2, "multiBigwigSummary/", x), recursive = T)
    dir.create(paste0(path2, "multiBigwigSummary/", x, "/zipped"), recursive = T)
    dir.create(paste0(path2, "multiBigwigSummary/", x, "/tabular"), recursive = T)
  })
}

# Run deepTools commands
print(paste0("Summarizing ", names(reg_type)[i], " regions for ",
             names(bw_files)[j], " (", j, " of ", length(bw_files), ")"))

system(paste0("multiBigwigSummary",
              " BED-file",
              " --bwfiles ", paste(bw_files[[j]], collapse = " "),
              " --BED ", reg_type[[i]],
              " --outFileName ", paste0(path2, "multiBigwigSummary/", names(bw_files)[j], "/zipped/", names(reg_type)[i] ,"_summarized_scores.npz"),
              " --outRawCounts ", paste0(path2, "multiBigwigSummary/", names(bw_files)[j], "/tabular/", names(reg_type)[i],"_summarized_scores.tab")
), wait = T)

system(paste0("multiBigwigSummary",
             " bins",
             " --bwfiles ", paste(bw_files[[j]], collapse = " "),
             " --outFileName ", paste0(path2, "multiBigwigSummary/", names(bw_files)[j], "/zipped/genomewide_summarized_scores.npz"),
             " --outRawCounts ", paste0(path2, "multiBigwigSummary/", names(bw_files)[j], "/tabular/genomewide_summarized_scores.tab"),
             " --numberOfProcessors 'max'"
), wait = T)

system(paste0("plotCorrelation",
             " -in ", paste0(path2, "multiBigwigSummary/", names(bw_files)[j], "/zipped/genomewide_summarized_scores.npz"),
             " --corMethod spearman --skipZeros --removeOutliers",
             " --plotHeight 15 --plotWidth 15",
             " --whatToPlot heatmap --colorMap RdYlBu --plotNumbers",
             " -o ", paste0(path2, "multiBigwigSummary/", names(bw_files)[j], "/zipped/genomewide_summarized_SpearmanCor.pdf")
), wait = T)

system(paste0("plotCorrelation",
             " -in ", paste0(path2, "multiBigwigSummary/", names(bw_files)[j], "/zipped/genomewide_summarized_scores.npz"),
             " --corMethod spearman --skipZeros --removeOutliers",
             " --plotHeight 15 --plotWidth 15",
             " --whatToPlot scatterplot",
             " -o ", paste0(path2, "multiBigwigSummary/", names(bw_files)[j], "/zipped/genomewide_summarized_SpearmanCor_Scatter.png")
), wait = T)
