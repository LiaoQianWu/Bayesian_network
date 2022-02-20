library("igraph")
library("visNetwork")
library("threejs")
library("htmlwidgets")

path <- "/path/to/dir/"
data_path <- "/path/to/summarizedMatrix/"

# Create output directory
if(!dir.exists(paste0(path, "interactiveViz"))){
  dir.create(paste0(path, "interactiveViz"), recursive = T)
}

# Read BN object
bn_network <- readRDS(paste0(path, "bayesianNetworksDiscretized_MAP/hct116_cpgrich_MU.RDS"))

# Read summarized matrix to calculate correlation between nodes
dat <- readRDS(paste0(data_path, "All_Consolidated_HCT116/hct116_cpgrich_MU_allEpigenomes.RDS"))
dat <- dat$Replicate1

# Compute correlation between nodes
bn_arcs <- bn_network$arcs
cor_nodes <- apply(bn_arcs, 1, function(x){
 a <- dat[, which(colnames(dat) == x[1])]
 b <- dat[, which(colnames(dat) == x[2])]
 c <- cor(a, b, method="spearman")
 return(c)})

# Create data frame
bn_nodes <- data.frame(id = names(bn_network$nodes))
bn_edges <- data.frame(from = bn_network$arc_strength$from, to = bn_network$arc_strength$to, weight = cor_nodes)

# Assign colors to each link
col_arcs <- rep("forestgreen", length(cor_nodes))
col_arcs[which(cor_nodes < 0)] <- "firebrick"

# Conduct visNetwork to generate interactive network
bn_nodes$shape <- "circle"
bn_nodes$label <- bn_nodes$id
bn_nodes$borderWidth <- 2
bn_nodes$color.background <- "aliceblue"
bn_nodes$color.border <- "white"
bn_nodes$color.highlight.background <- "blanchedalmond"
bn_nodes$color.highlight.border <- "white"
bn_nodes$shadow <- T

bn_edges$width <- abs(bn_edges$weight) * 10
bn_edges$color <- col_arcs
bn_edges$arrows <- "to"

net_vis <- visNetwork(bn_nodes, bn_edges, width = "100%", height = "700px") %>%
  visOptions(highlightNearest = T, selectedBy = "label") %>%
  visLayout(randomSeed = 123, improvedLayout = T)

# Save network as html file
visSave(net_vis, file = paste0(path, "interactiveViz/hct116_cpgrich_MU.html"))
