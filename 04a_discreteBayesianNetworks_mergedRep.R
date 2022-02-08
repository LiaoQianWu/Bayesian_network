library(bnlearn)
library(parallel)

path <- "/path/to/dir/"

# Define region of interest
regType <- c(hct116_cpgrich_MU = "hct116_cpgrich_MU", hct116_cpgpoor_MU = "hct116_cpgpoor_MU")
intervention <- c(dko1_cpgrich_MU = "dko1_cpgrich_MU", dko1_cpgpoor_MU = "dko1_cpgpoor_MU")

################################################################################
# MAP Bayesian
################################################################################

# Create output directory
if (!dir.exists(paste0(path, "bayesianNetworksDiscretized_MAP"))){
  dir.create(paste0(path, "bayesianNetworksDiscretized_MAP"), recursive = T)}

################################################################################
# Use mixed observational and interventional data to learn Bayesian networks
################################################################################

for(i in 1:length(regType)){
  
  # Load observational data (discretized consolidated matrices)
  dat <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/", regType[[i]], "_discretized.RDS"))
  
  # Replicate 1
  dat_r1 <- dat[["Replicate1"]]
  dat_r1 <- apply(dat_r1, 2, as.character)
  dat_r1 <- data.frame(dat_r1, stringsAsFactors = T)
  dat_r1$INT <- as.factor(0)
  
  # Replicate 2
  dat_r2 <- dat[["Replicate2"]]
  dat_r2 <- apply(dat_r2, 2, as.character)
  dat_r2 <- data.frame(dat_r2, stringsAsFactors = T)
  dat_r2$INT <- as.factor(0)
  
  # Load interventional data (discretized consolidated matrices)
  int_dat <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/", intervention[[i]], "_discretized.RDS"))
  
  # Replicate 1
  int_dat_r1 <- int_dat[["Replicate1"]]
  # Manually set methylation column to the lowest bin i.e. 1.
  int_dat_r1[,2] <- rep(1, nrow(int_dat_r1))
  int_dat_r1 <- apply(int_dat_r1, 2, as.character)
  int_dat_r1 <- data.frame(int_dat_r1, stringsAsFactors = T)
  int_dat_r1$INT <- as.factor(2)
  
  # Replicate 2
  int_dat_r2 <- int_dat[["Replicate2"]]
  # Manually set methylation column to the lowest bin i.e. 1.
  int_dat_r2[,2] <- rep(1, nrow(int_dat_r2))
  int_dat_r2 <- apply(int_dat_r2, 2, as.character)
  int_dat_r2 <- data.frame(int_dat_r2, stringsAsFactors = T)
  int_dat_r2$INT <- as.factor(2)
  
  # Combine observational and interventional data and create interventional index list for latter use
  r1 <- rbind(dat_r1, int_dat_r1)
  r2 <- rbind(dat_r2, int_dat_r2)
  INT <- sapply(1:10, function(x){which(r1$INT == x)})
  r1 <- r1[,-which(colnames(r1) == "INT")]
  r2 <- r2[,-which(colnames(r2) == "INT")]
  names(INT) <- names(r1)
  
  # Define blacklist to exclude direction from expression to any other variables
  bl <- data.frame(from = "Expression", to = colnames(r1))
  bl <- bl[-which(bl$to == "Expression"), ]
  
  # Structure learning
  # Start from a set of random graphs
  nodes <- names(r1)
  init <- random.graph(nodes = nodes, method = 'melancon', num = 500, burn.in = 10^5, every = 100)
  
  # Tune each random graphs with data
  netlist_r1 <- lapply(init, function(net){tabu(r1, score = 'mbde', exp = INT, start = net, tabu = 50, blacklist = bl)})
  netlist_r2 <- lapply(init, function(net){tabu(r2, score = 'mbde', exp = INT, start = net, tabu = 50, blacklist = bl)})
  
  # Take resulting list to compute arc strengths
  strength_r1 <- custom.strength(netlist_r1, nodes = nodes, cpdag = F)
  strength_r2 <- custom.strength(netlist_r2, nodes = nodes, cpdag = F)
  strength <- mean(strength_r1, strength_r2)
  
  # Model averaging
  strength_filt <- strength[strength$strength > 0.85 & strength$direction > 0.5,]
  dag_average <- averaged.network(strength_filt, threshold = 0.85)
  dag_average$learning$blacklist <- bl
  #dag_average$arc_strength <- strength_filt
  dag_average$arc_strength <- strength #More flexible for downstream operations
  
  # Save results of structure learning
  saveRDS(dag_average, paste0(path, "bayesianNetworksDiscretized_MAP/", regType[[i]], ".RDS"))
}

################################################################################
# Frequentist (bootstrap)
################################################################################

cores_cnt = makeCluster(60)

# Create output directory
if (!dir.exists(paste0(path, "bayesianNetworksDiscretized_bootstrap"))){
  dir.create(paste0(path, "bayesianNetworksDiscretized_bootstrap"), recursive = T)}

################################################################################
# Use mixed observational and interventional data to learn Bayesian networks
################################################################################

for(i in 1:length(regType)){
  
  # Load observational data (discretized consolidated matrices)
  dat <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/", regType[[i]], "_discretized.RDS"))
  
  # Replicate 1
  dat_r1 <- dat[["Replicate1"]]
  dat_r1 <- apply(dat_r1, 2, as.character)
  dat_r1 <- data.frame(dat_r1, stringsAsFactors = T)
  dat_r1$INT <- as.factor(0)
  
  # Replicate 2
  dat_r2 <- dat[["Replicate2"]]
  dat_r2 <- apply(dat_r2, 2, as.character)
  dat_r2 <- data.frame(dat_r2, stringsAsFactors = T)
  dat_r2$INT <- as.factor(0)
  
  # Load interventional data (discretized consolidated matrices)
  int_dat <- readRDS(paste0(path, "summarizedMatrix/All_Consolidated_Epigenomes_Discretized/", intervention[[i]], "_discretized.RDS"))
  
  # Replicate 1
  int_dat_r1 <- int_dat[["Replicate1"]]
  # Manually set methylation column to the lowest bin i.e. 1.
  int_dat_r1[,2] <- rep(1, nrow(int_dat_r1))
  int_dat_r1 <- apply(int_dat_r1, 2, as.character)
  int_dat_r1 <- data.frame(int_dat_r1, stringsAsFactors = T)
  int_dat_r1$INT <- as.factor(2)
  
  # Replicate 2
  int_dat_r2 <- int_dat[["Replicate2"]]
  # Manually set methylation column to the lowest bin i.e. 1.
  int_dat_r2[,2] <- rep(1, nrow(int_dat_r2))
  int_dat_r2 <- apply(int_dat_r2, 2, as.character)
  int_dat_r2 <- data.frame(int_dat_r2, stringsAsFactors = T)
  int_dat_r2$INT <- as.factor(2)
  
  # Combine observational and interventional data and create interventional index list for latter use
  r1 <- rbind(dat_r1, int_dat_r1)
  r2 <- rbind(dat_r2, int_dat_r2)
  INT <- sapply(1:10, function(x){which(r1$INT == x)})
  r1 <- r1[,-which(colnames(r1) == "INT")]
  r2 <- r2[,-which(colnames(r2) == "INT")]
  names(INT) <- names(r1)
  
  # Define blacklist to exclude direction from expression to any other variables
  bl <- data.frame(from = "Expression", to = colnames(r1))
  bl <- bl[-which(bl$to == "Expression"), ]
  
  # Structure learning
  strength_r1 = boot.strength(r1, algorithm = "tabu",
                              algorithm.args = list(score = "mbde", exp = INT, tabu = 50, blacklist = bl),
                              cluster = cores_cnt, R = 500, cpdag = F, m = round(nrow(dat_r1) * 0.85))
  strength_r2 = boot.strength(r2, algorithm = "tabu",
                              algorithm.args = list(score = "mbde", exp = INT, tabu = 50, blacklist = bl), 
                              cluster = cores_cnt, R = 500, cpdag = F, m = round(nrow(dat_r1) * 0.85))
  strength <- mean(strength_r1, strength_r2)
  
  # Model averaging
  strength_filt <- strength[strength$strength > 0.90 & strength$direction > 0.5,]
  dag_average <- averaged.network(strength, threshold = 0.85)
  dag_average$learning$blacklist <- bl
  dag_average$arc_strength <- strength
  
  # Save results of structure learning
  saveRDS(dag_average, paste0(path, "bayesianNetworksDiscretized_bootstrap/", regType[[i]], ".RDS"))
}
stopCluster(cores_cnt)
