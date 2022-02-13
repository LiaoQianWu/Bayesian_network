library(bnlearn)
library(igraph)
library(kableExtra)

# Import igraph2 for acommodating multiple arrow size settings
# Code was taken from jevansbio (https://github.com/jevansbio/igraphhack)
source("/path/to/igraphhack-master/igraphplot2.R")
environment(plot.igraph2) <- asNamespace("igraph")
environment(igraph.Arrows2) <- asNamespace("igraph")

path <- "/path/to/dir/"

################################################################################
# WITH INTERVENTIONAL DATA
# Bayesian network info
################################################################################

# Access to bn objects
hct116_cpgrich_MU <- readRDS(paste0(path, "bayesianNetworksDiscretized_bootstrap/hct116_cpgrich_MU.RDS"))
hct116_cpgpoor_MU <- readRDS(paste0(path, "bayesianNetworksDiscretized_bootstrap/hct116_cpgpoor_MU.RDS"))

################################################################################
# Continuous matrices
################################################################################

# Access to prior feature information for plot settings
dat_hct116_cpgrich_MU <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_HCT116/hct116_cpgrich_MU_allEpigenomes.RDS"))
dat_hct116_cpgrich_MU_r1 <- as.data.frame(dat_hct116_cpgrich_MU$Replicate1)
dat_hct116_cpgrich_MU_r2 <- as.data.frame(dat_hct116_cpgrich_MU$Replicate2)

dat_dko1_cpgrich_MU <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_DKO1/dko1_cpgrich_MU_allEpigenomes.RDS"))
dat_dko1_cpgrich_MU_r1 <- as.data.frame(dat_dko1_cpgrich_MU$Replicate1)
dat_dko1_cpgrich_MU_r2 <- as.data.frame(dat_dko1_cpgrich_MU$Replicate2)

dat_hct116_cpgpoor_MU <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_HCT116/hct116_cpgpoor_MU_allEpigenomes.RDS"))
dat_hct116_cpgpoor_MU_r1 <- as.data.frame(dat_hct116_cpgpoor_MU$Replicate1)
dat_hct116_cpgpoor_MU_r2 <- as.data.frame(dat_hct116_cpgpoor_MU$Replicate2)

dat_dko1_cpgpoor_MU <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_DKO1/dko1_cpgpoor_MU_allEpigenomes.RDS"))
dat_dko1_cpgpoor_MU_r1 <- as.data.frame(dat_dko1_cpgpoor_MU$Replicate1)
dat_dko1_cpgpoor_MU_r2 <- as.data.frame(dat_dko1_cpgpoor_MU$Replicate2)

################################################################################
# Correlation between arcs
################################################################################

########################
# HCT116 CGI promoter MU
########################
hct116_cpgrich_MU_arcs <- as.data.frame(hct116_cpgrich_MU$arcs)

# Investigate relationship between HCT116 CGI promoter MU rep1&rep2 and DKO1 CGI promoter MU rep1&rep2
# because we used DKO1 as interventional data to improve structure learning
cor_hct116_cpgrich_MU_r1 <- apply(hct116_cpgrich_MU_arcs, 1, function(x){
  a <- dat_hct116_cpgrich_MU_r1[, which(colnames(dat_hct116_cpgrich_MU_r1) == x[1])]
  b <- dat_hct116_cpgrich_MU_r1[, which(colnames(dat_hct116_cpgrich_MU_r1) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

cor_hct116_cpgrich_MU_r2 <- apply(hct116_cpgrich_MU_arcs, 1, function(x){
  a <- dat_hct116_cpgrich_MU_r2[, which(colnames(dat_hct116_cpgrich_MU_r2) == x[1])]
  b <- dat_hct116_cpgrich_MU_r2[, which(colnames(dat_hct116_cpgrich_MU_r2) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

cor_dko1_cpgrich_MU_r1 <- apply(hct116_cpgrich_MU_arcs, 1, function(x){
  a <- dat_dko1_cpgrich_MU_r1[, which(colnames(dat_dko1_cpgrich_MU_r1) == x[1])]
  b <- dat_dko1_cpgrich_MU_r1[, which(colnames(dat_dko1_cpgrich_MU_r1) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

cor_dko1_cpgrich_MU_r2 <- apply(hct116_cpgrich_MU_arcs, 1, function(x){
  a <- dat_dko1_cpgrich_MU_r2[, which(colnames(dat_dko1_cpgrich_MU_r2) == x[1])]
  b <- dat_dko1_cpgrich_MU_r2[, which(colnames(dat_dko1_cpgrich_MU_r2) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

# Level correlation >, <, and = 0
cor_a <- vector()
for(i in 1:length(cor_hct116_cpgrich_MU_r1)){
  if(cor_hct116_cpgrich_MU_r1[i] > 0){
    cor_a[i] <- 1
  }
  if(cor_hct116_cpgrich_MU_r1[i] < 0){
    cor_a[i] <- -1
  }
  if(cor_hct116_cpgrich_MU_r1[i] == 0){
    cor_a[i] <- 0
  }
}

cor_b <- vector()
for(i in 1:length(cor_hct116_cpgrich_MU_r2)){
  if(cor_hct116_cpgrich_MU_r2[i] > 0){
    cor_b[i] <- 1
  }
  if(cor_hct116_cpgrich_MU_r2[i] < 0){
    cor_b[i] <- -1
  }
  if(cor_hct116_cpgrich_MU_r2[i] == 0){
    cor_b[i] <- 0
  }
}

cor_c <- vector()
for(i in 1:length(cor_dko1_cpgrich_MU_r1)){
  if(cor_dko1_cpgrich_MU_r1[i] > 0){
    cor_c[i] <- 1
  }
  if(cor_dko1_cpgrich_MU_r1[i] < 0){
    cor_c[i] <- -1
  }
  if(cor_dko1_cpgrich_MU_r1[i] == 0){
    cor_c[i] <- 0
  }
}

cor_d <- vector()
for(i in 1:length(cor_dko1_cpgrich_MU_r2)){
  if(cor_dko1_cpgrich_MU_r2[i] > 0){
    cor_d[i] <- 1
  }
  if(cor_dko1_cpgrich_MU_r2[i] < 0){
    cor_d[i] <- -1
  }
  if(cor_dko1_cpgrich_MU_r2[i] == 0){
    cor_d[i] <- 0
  }
}

# Append levels of correlation to arc matrix
hct116_cpgrich_MU_arcs$hct116_r1 <- cor_a
hct116_cpgrich_MU_arcs$hct116_r2 <- cor_b
hct116_cpgrich_MU_arcs$dko1_r1 <- cor_c
hct116_cpgrich_MU_arcs$dko1_r2 <- cor_d

rm(cor_hct116_cpgrich_MU_r1, cor_hct116_cpgrich_MU_r2,
   cor_dko1_cpgrich_MU_r1, cor_dko1_cpgrich_MU_r2,
   cor_a, cor_b, cor_c, cor_d)

############################
# HCT116 non-CGI promoter MU
############################
hct116_cpgpoor_MU_arcs <- as.data.frame(hct116_cpgpoor_MU$arcs)

# Investigate relationship between HCT116 non-CGI promoter MU rep1&rep2 and DKO1 non-CGI promoter MU rep1&rep2
# because we used DKO1 as interventional data to improve structure learning
cor_hct116_cpgpoor_MU_r1 <- apply(hct116_cpgpoor_MU_arcs, 1, function(x){
  a <- dat_hct116_cpgpoor_MU_r1[, which(colnames(dat_hct116_cpgpoor_MU_r1) == x[1])]
  b <- dat_hct116_cpgpoor_MU_r1[, which(colnames(dat_hct116_cpgpoor_MU_r1) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

cor_hct116_cpgpoor_MU_r2 <- apply(hct116_cpgpoor_MU_arcs, 1, function(x){
  a <- dat_hct116_cpgpoor_MU_r2[, which(colnames(dat_hct116_cpgpoor_MU_r2) == x[1])]
  b <- dat_hct116_cpgpoor_MU_r2[, which(colnames(dat_hct116_cpgpoor_MU_r2) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

cor_dko1_cpgpoor_MU_r1 <- apply(hct116_cpgpoor_MU_arcs, 1, function(x){
  a <- dat_dko1_cpgpoor_MU_r1[, which(colnames(dat_dko1_cpgpoor_MU_r1) == x[1])]
  b <- dat_dko1_cpgpoor_MU_r1[, which(colnames(dat_dko1_cpgpoor_MU_r1) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

cor_dko1_cpgpoor_MU_r2 <- apply(hct116_cpgpoor_MU_arcs, 1, function(x){
  a <- dat_dko1_cpgpoor_MU_r2[, which(colnames(dat_dko1_cpgpoor_MU_r2) == x[1])]
  b <- dat_dko1_cpgpoor_MU_r2[, which(colnames(dat_dko1_cpgpoor_MU_r2) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

# Level correlation >, <, and = 0
cor_a <- vector()
for(i in 1:length(cor_hct116_cpgpoor_MU_r1)){
  if(cor_hct116_cpgpoor_MU_r1[i] > 0){
    cor_a[i] <- 1
  }
  if(cor_hct116_cpgpoor_MU_r1[i] < 0){
    cor_a[i] <- -1
  }
  if(cor_hct116_cpgpoor_MU_r1[i] == 0){
    cor_a[i] <- 0
  }
}

cor_b <- vector()
for(i in 1:length(cor_hct116_cpgpoor_MU_r2)){
  if(cor_hct116_cpgpoor_MU_r2[i] > 0){
    cor_b[i] <- 1
  }
  if(cor_hct116_cpgpoor_MU_r2[i] < 0){
    cor_b[i] <- -1
  }
  if(cor_hct116_cpgpoor_MU_r2[i] == 0){
    cor_b[i] <- 0
  }
}

cor_c <- vector()
for(i in 1:length(cor_dko1_cpgpoor_MU_r1)){
  if(cor_dko1_cpgpoor_MU_r1[i] > 0){
    cor_c[i] <- 1
  }
  if(cor_dko1_cpgpoor_MU_r1[i] < 0){
    cor_c[i] <- -1
  }
  if(cor_dko1_cpgpoor_MU_r1[i] == 0){
    cor_c[i] <- 0
  }
}

cor_d <- vector()
for(i in 1:length(cor_dko1_cpgpoor_MU_r2)){
  if(cor_dko1_cpgpoor_MU_r2[i] > 0){
    cor_d[i] <- 1
  }
  if(cor_dko1_cpgpoor_MU_r2[i] < 0){
    cor_d[i] <- -1
  }
  if(cor_dko1_cpgpoor_MU_r2[i] == 0){
    cor_d[i] <- 0
  }
}

# Append levels of correlation to arc matrix
hct116_cpgpoor_MU_arcs$hct116_r1 <- cor_a
hct116_cpgpoor_MU_arcs$hct116_r2 <- cor_b
hct116_cpgpoor_MU_arcs$dko1_r1 <- cor_c
hct116_cpgpoor_MU_arcs$dko1_r2 <- cor_d

rm(cor_hct116_cpgpoor_MU_r1, cor_hct116_cpgpoor_MU_r2,
   cor_dko1_cpgpoor_MU_r1, cor_dko1_cpgpoor_MU_r2,
   cor_a, cor_b, cor_c, cor_d)

################################################################################
# plot network structures
################################################################################

regType <- list(hct116_cpgrich_MU, hct116_cpgpoor_MU)

titleName <- c("HCT116 CGI Promoter MU", "HCT116 Non-CGI Promoter MU")

# Create two-column matrix for highlight argument in strength.plot
neg_cor_hct116_cpgrich_MU <- hct116_cpgrich_MU_arcs[which(hct116_cpgrich_MU_arcs$hct116_r1 == -1), c(1,2)]
neg_cor_hct116_cpgpoor_MU <- hct116_cpgpoor_MU_arcs[which(hct116_cpgpoor_MU_arcs$hct116_r1 == -1), c(1,2)]
hl_arcs <- list(neg_cor_hct116_cpgrich_MU, neg_cor_hct116_cpgpoor_MU)

for(i in 1:length(regType)){
  bn <- regType[[i]]
  strength.plot(x = bn,
                strength = bn$arc_strength[bn$arc_strength$strength > 0.85 & bn$arc_strength$direction > 0.5,],
                main = titleName[i],
                highlight = list(nodes = c("Methylation", nbr(x=bn, node = "Methylation")), arcs = hl_arcs[[i]], lty = "longdash", fill = "pink"))
}

################################################################################
# igraph
################################################################################

bn_object_list <- list(hct116_cpgrich_MU,
                       hct116_cpgpoor_MU)
dat_list <- list(dat_hct116_cpgrich_MU_r1,
                 dat_hct116_cpgpoor_MU_r1)

title_list <- c("HCT116 CGI Promoter MU", "HCT116 Non-CGI Promoter MU")

for(i in 1:length(bn_object_list)){
  arc_strength <- bn_object_list[[i]]$arc_strength
  bn_arcs <- as.matrix(arc_strength[arc_strength$strength > 0.85 & arc_strength$direction != 0, 1:2])
  bn_arc_direction <- arc_strength[arc_strength$strength > 0.85 & arc_strength$direction != 0, 4]
  
  # Investigate into direction probabilities
  print(title_list[i])
  print(arc_strength[arc_strength$strength > 0.85 & arc_strength$direction != 0, c(1,2,4)])
  
  dat <- dat_list[[i]]
  cor_nodes <- apply(bn_arcs, 1, function(x){
    a <- dat[, which(colnames(dat) == x[1])]
    b <- dat[, which(colnames(dat) == x[2])]
    c <- cor(a, b, method="spearman")
    return(c)})
  
  col_arcs <- rep("forestgreen", length(cor_nodes))
  col_arcs[which(cor_nodes < 0)] <- "firebrick"
  
  bn_network <- graph_from_edgelist(bn_arcs, directed = T)
  E(bn_network)$rho <- cor_nodes
  E(bn_network)$direction <- bn_arc_direction
  
  plot.igraph2(bn_network, vertex.shape="none", vertex.label.font=2, vertex.label.cex=0.8,
               edge.color=col_arcs, edge.width=0.5+abs(E(bn_network)$rho)*10,
               edge.arrow.size=0.5+E(bn_network)$direction*1.8, layout=layout_as_star, main=title_list[i])
  legend(-1.65, -1.15, legend=c("Positive", "Negative"), col=c("forestgreen", "firebrick"),
         lty=1, lwd=2, cex = 0.8, title="Correlation", bty="n")
  legend(-0.95, -1.25, legend=c(0.1, 0.2, 0.3, 0.4), col="black",
         lty=1, lwd=0.5+c(0.1, 0.2, 0.3, 0.4)*10, cex=0.8, bty="n")
  legend(-0.47, -1.25, legend=c(0.5, 0.6), col="black",
         lty=1, lwd=0.5+c(0.5, 0.6)*10, cex=0.8, bty="n")
  legend(0.25, -1.15, legend=c(0.2, 0.4, 0.6), col="black", title="Direction probability",
         pch=17, pt.cex=0.5+c(0.2, 0.4, 0.6)*1.8, cex=0.8, bty="n")
  legend(0.9, -1.265, legend=c(0.8), col="black",
         pch=17, pt.cex=0.5+c(0.8)*1.8, cex=0.8, bty="n")
  
  # igraph
  # plot(bn_network, vertex.shape="none", vertex.label.font=2, vertex.label.cex=0.8,
  #      edge.color=col_arcs, edge.width=0.5+abs(E(bn_network)$rho)*10,
  #      edge.arrow.size=1.5, layout=layout_as_star, main=title_list[i])
}


################################################################################
# WITHOUT INTERVENTIONAL DATA
# Bayesian network info
################################################################################

# Access to bn objects
no_int_hct116_cpgrich_MU <- readRDS(paste0(path, "bayesianNetworksDiscretized_MAP/no_int_hct116_cpgrich_MU.RDS"))
no_int_hct116_cpgpoor_MU <- readRDS(paste0(path, "bayesianNetworksDiscretized_MAP/no_int_hct116_cpgpoor_MU.RDS"))

################################################################################
# Continuous matrices
################################################################################

# Access to prior feature information for plot settings
dat_hct116_cpgrich_MU <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_HCT116/hct116_cpgrich_MU_allEpigenomes.RDS"))
dat_hct116_cpgrich_MU_r1 <- as.data.frame(dat_hct116_cpgrich_MU$Replicate1)
dat_hct116_cpgrich_MU_r2 <- as.data.frame(dat_hct116_cpgrich_MU$Replicate2)

dat_hct116_cpgpoor_MU <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_HCT116/hct116_cpgpoor_MU_allEpigenomes.RDS"))
dat_hct116_cpgpoor_MU_r1 <- as.data.frame(dat_hct116_cpgpoor_MU$Replicate1)
dat_hct116_cpgpoor_MU_r2 <- as.data.frame(dat_hct116_cpgpoor_MU$Replicate2)

################################################################################
# Correlation between arcs
################################################################################

########################
# HCT116 CGI promoter MU
########################
hct116_cpgrich_MU_arcs <- as.data.frame(no_int_hct116_cpgrich_MU$arcs)

# Investigate relationship between HCT116 CGI promoter MU rep1&rep2
cor_hct116_cpgrich_MU_r1 <- apply(hct116_cpgrich_MU_arcs, 1, function(x){
  a <- dat_hct116_cpgrich_MU_r1[, which(colnames(dat_hct116_cpgrich_MU_r1) == x[1])]
  b <- dat_hct116_cpgrich_MU_r1[, which(colnames(dat_hct116_cpgrich_MU_r1) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

cor_hct116_cpgrich_MU_r2 <- apply(hct116_cpgrich_MU_arcs, 1, function(x){
  a <- dat_hct116_cpgrich_MU_r2[, which(colnames(dat_hct116_cpgrich_MU_r2) == x[1])]
  b <- dat_hct116_cpgrich_MU_r2[, which(colnames(dat_hct116_cpgrich_MU_r2) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

# Level correlation >, <, and = 0
cor_a <- vector()
for(i in 1:length(cor_hct116_cpgrich_MU_r1)){
  if(cor_hct116_cpgrich_MU_r1[i] > 0){
    cor_a[i] <- 1
  }
  if(cor_hct116_cpgrich_MU_r1[i] < 0){
    cor_a[i] <- -1
  }
  if(cor_hct116_cpgrich_MU_r1[i] == 0){
    cor_a[i] <- 0
  }
}

cor_b <- vector()
for(i in 1:length(cor_hct116_cpgrich_MU_r2)){
  if(cor_hct116_cpgrich_MU_r2[i] > 0){
    cor_b[i] <- 1
  }
  if(cor_hct116_cpgrich_MU_r2[i] < 0){
    cor_b[i] <- -1
  }
  if(cor_hct116_cpgrich_MU_r2[i] == 0){
    cor_b[i] <- 0
  }
}

# Append levels of correlation to arc matrix
hct116_cpgrich_MU_arcs$hct116_r1 <- cor_a
hct116_cpgrich_MU_arcs$hct116_r2 <- cor_b

rm(cor_hct116_cpgrich_MU_r1, cor_hct116_cpgrich_MU_r2, cor_a, cor_b)

############################
# HCT116 non-CGI promoter MU
############################
hct116_cpgpoor_MU_arcs <- as.data.frame(no_int_hct116_cpgpoor_MU$arcs)

# Investigate relationship between HCT116 non-CGI promoter MU rep1&rep2
cor_hct116_cpgpoor_MU_r1 <- apply(hct116_cpgpoor_MU_arcs, 1, function(x){
  a <- dat_hct116_cpgpoor_MU_r1[, which(colnames(dat_hct116_cpgpoor_MU_r1) == x[1])]
  b <- dat_hct116_cpgpoor_MU_r1[, which(colnames(dat_hct116_cpgpoor_MU_r1) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

cor_hct116_cpgpoor_MU_r2 <- apply(hct116_cpgpoor_MU_arcs, 1, function(x){
  a <- dat_hct116_cpgpoor_MU_r2[, which(colnames(dat_hct116_cpgpoor_MU_r2) == x[1])]
  b <- dat_hct116_cpgpoor_MU_r2[, which(colnames(dat_hct116_cpgpoor_MU_r2) == x[2])]
  c <- cor(a, b, method="spearman")
  return(c)})

# Level correlation >, <, and = 0
cor_a <- vector()
for(i in 1:length(cor_hct116_cpgpoor_MU_r1)){
  if(cor_hct116_cpgpoor_MU_r1[i] > 0){
    cor_a[i] <- 1
  }
  if(cor_hct116_cpgpoor_MU_r1[i] < 0){
    cor_a[i] <- -1
  }
  if(cor_hct116_cpgpoor_MU_r1[i] == 0){
    cor_a[i] <- 0
  }
}

cor_b <- vector()
for(i in 1:length(cor_hct116_cpgpoor_MU_r2)){
  if(cor_hct116_cpgpoor_MU_r2[i] > 0){
    cor_b[i] <- 1
  }
  if(cor_hct116_cpgpoor_MU_r2[i] < 0){
    cor_b[i] <- -1
  }
  if(cor_hct116_cpgpoor_MU_r2[i] == 0){
    cor_b[i] <- 0
  }
}

# Append levels of correlation to arc matrix
hct116_cpgpoor_MU_arcs$hct116_r1 <- cor_a
hct116_cpgpoor_MU_arcs$hct116_r2 <- cor_b

rm(cor_hct116_cpgpoor_MU_r1, cor_hct116_cpgpoor_MU_r2, cor_a, cor_b)

################################################################################
# plot network structures
################################################################################

regType <- list(no_int_hct116_cpgrich_MU, no_int_hct116_cpgpoor_MU)

titleName <- c("HCT116 CGI Promoter MU (without interventional data)", "HCT116 Non-CGI Promoter MU (without interventional data)")

# Create two-column matrix for highlight argument in strength.plot
neg_cor_hct116_cpgrich_MU <- hct116_cpgrich_MU_arcs[which(hct116_cpgrich_MU_arcs$hct116_r1 == -1), c(1,2)]
neg_cor_hct116_cpgpoor_MU <- hct116_cpgpoor_MU_arcs[which(hct116_cpgpoor_MU_arcs$hct116_r1 == -1), c(1,2)]
hl_arcs <- list(neg_cor_hct116_cpgrich_MU, neg_cor_hct116_cpgpoor_MU)

for(i in 1:length(regType)){
  bn <- regType[[i]]
  strength.plot(x = bn,
                strength = bn$arc_strength[bn$arc_strength$strength > 0.75 & bn$arc_strength$direction > 0.5,],
                main = titleName[i],
                highlight = list(nodes = c("Methylation", nbr(x=bn, node = "Methylation")), arcs = hl_arcs[[i]], lty = "longdash", fill = "pink"))
}

################################################################################
# igraph
################################################################################

bn_object_list <- list(no_int_hct116_cpgrich_MU,
                       no_int_hct116_cpgpoor_MU)
dat_list <- list(dat_hct116_cpgrich_MU_r1,
                 dat_hct116_cpgpoor_MU_r1)

title_list <- c("HCT116 CGI Promoter MU (without interventional data)", "HCT116 Non-CGI Promoter MU (without interventional data)")

for(i in 1:length(bn_object_list)){
  arc_strength <- bn_object_list[[i]]$arc_strength
  bn_arcs <- as.matrix(arc_strength[arc_strength$strength > 0.85 & arc_strength$direction != 0, 1:2])
  bn_arc_direction <- arc_strength[arc_strength$strength > 0.85 & arc_strength$direction != 0, 4]
  
  # Investigate into direction probabilities
  print(title_list[i])
  print(arc_strength[arc_strength$strength > 0.85 & arc_strength$direction != 0, c(1,2,4)])
  
  dat <- dat_list[[i]]
  cor_nodes <- apply(bn_arcs, 1, function(x){
    a <- dat[, which(colnames(dat) == x[1])]
    b <- dat[, which(colnames(dat) == x[2])]
    c <- cor(a, b, method="spearman")
    return(c)})
  
  col_arcs <- rep("forestgreen", length(cor_nodes))
  col_arcs[which(cor_nodes < 0)] <- "firebrick"
  
  bn_network <- graph_from_edgelist(bn_arcs, directed = T)
  E(bn_network)$rho <- cor_nodes
  E(bn_network)$direction <- bn_arc_direction
  
  plot.igraph2(bn_network, vertex.shape="none", vertex.label.font=2, vertex.label.cex=0.8,
               edge.color=col_arcs, edge.width=0.5+abs(E(bn_network)$rho)*10,
               edge.arrow.size=0.5+E(bn_network)$direction*1.8, layout=layout_as_star, main=title_list[i])
  legend(-1.65, -1.15, legend=c("Positive", "Negative"), col=c("forestgreen", "firebrick"),
         lty=1, lwd=2, cex = 0.8, title="Correlation", bty="n")
  legend(-0.95, -1.25, legend=c(0.1, 0.2, 0.3, 0.4), col="black",
         lty=1, lwd=0.5+c(0.1, 0.2, 0.3, 0.4)*10, cex=0.8, bty="n")
  legend(-0.47, -1.25, legend=c(0.5, 0.6), col="black",
         lty=1, lwd=0.5+c(0.5, 0.6)*10, cex=0.8, bty="n")
  legend(0.25, -1.15, legend=c(0.2, 0.4, 0.6), col="black", title="Direction probability",
         pch=17, pt.cex=0.5+c(0.2, 0.4, 0.6)*1.8, cex=0.8, bty="n")
  legend(0.9, -1.265, legend=c(0.8), col="black",
         pch=17, pt.cex=0.5+c(0.8)*1.8, cex=0.8, bty="n")
}
