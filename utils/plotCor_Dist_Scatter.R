library(pheatmap)
library(reshape2)
library(ggplot2)
library(GGally)
library(RColorBrewer)

path <- "/path/to/summarizedMatrix/"

################################################################################
# Continuous data
################################################################################

hct116_cpgrich_MU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgrich_MU_allEpigenomes.RDS"))
hct116_cpgrich_UU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgrich_UU_allEpigenomes.RDS"))
hct116_cpgpoor_MU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgpoor_MU_allEpigenomes.RDS"))
hct116_cpgpoor_UU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgpoor_UU_allEpigenomes.RDS"))
dko1_cpgrich_MU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgrich_MU_allEpigenomes.RDS"))
dko1_cpgrich_UU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgrich_UU_allEpigenomes.RDS"))
dko1_cpgpoor_MU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgpoor_MU_allEpigenomes.RDS"))
dko1_cpgpoor_UU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgpoor_UU_allEpigenomes.RDS"))

r <- c("Replicate1", "Replicate2")
summarized_matrix_list <- list(hct116_cpgrich_MU, hct116_cpgrich_UU,
                               hct116_cpgpoor_MU, hct116_cpgpoor_UU,
                               dko1_cpgrich_MU, dko1_cpgrich_UU,
                               dko1_cpgpoor_MU, dko1_cpgpoor_UU)

title_list <- c("HCT116 CGI Promoter MU", "HCT116 CGI Promoter UU",
                "HCT116 Non-CGI Promoter MU", "HCT116 Non-CGI Promoter UU",
                "DKO1 CGI Promoter MU", "DKO1 CGI Promoter UU",
                "DKO1 Non-CGI Promoter MU", "DKO1 Non-CGI Promoter UU")

################################################################################
# Generate correlation heatmap using continuous data
################################################################################

# Create output directory
if(!dir.exists(paste0(path, "plotCorHeatmap"))){
  dir.create(paste0(path, "plotCorHeatmap"))
}

for(i in 1:length(summarized_matrix_list)){
  rep1 <- summarized_matrix_list[[i]]$Replicate1
  rep2 <- summarized_matrix_list[[i]]$Replicate2
  meth <- data.frame(Methylation = rep1[, 2])
  rep1 <- rep1[, -2]
  rep2 <- rep2[, -2]
  
  if(identical(colnames(rep1), colnames(rep2)) & identical(rownames(rep1), rownames(rep2))){
    colnames(rep1) <- paste0(colnames(rep1), "_Rep1")
    colnames(rep2) <- paste0(colnames(rep2), "_Rep2")
    if(sum(meth) != 0){
      combo <- cbind(meth, rep1, rep2) 
    }
    else{
      combo <- cbind(rep1, rep2)
    }
  }
  
  combo_cor <- cor(combo, method = "spearman")
  pheatmap(combo_cor, main = title_list[i], display_numbers = T, angle_col = "45")
}

################################################################################
# Display distribution of continuous data
################################################################################

color <- brewer.pal(10, "Set3")

# HCT116 CGI promoter MU
tmp_r1 <- melt(hct116_cpgrich_MU$Replicate1)
tmp_r2 <- melt(hct116_cpgrich_MU$Replicate2)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("Data Distribution of HCT116 CGI Promoter MU") +
  labs(x = "log2(value)") + guides(color = FALSE) +
  scale_fill_discrete(labels = c("Rep1", "Rep2")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

# HCT116 CGI promoter UU
tmp_r1 <- melt(hct116_cpgrich_UU$Replicate1)
tmp_r2 <- melt(hct116_cpgrich_UU$Replicate2)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("Data Distribution of HCT116 CGI Promoter UU") +
  labs(x = "log2(value)") + guides(color = FALSE) +
  scale_fill_discrete(labels = c("Rep1", "Rep2")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

# DKO1 CGI promoter MU
tmp_r1 <- melt(dko1_cpgrich_MU$Replicate1)
tmp_r2 <- melt(dko1_cpgrich_MU$Replicate2)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("Data Distribution of DKO1 CGI Promoter MU") +
  labs(x = "log2(value)") + guides(color = FALSE) +
  scale_fill_discrete(labels = c("Rep1", "Rep2")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

# DKO1 CGI promoter UU
tmp_r1 <- melt(dko1_cpgrich_UU$Replicate1)
tmp_r2 <- melt(dko1_cpgrich_UU$Replicate2)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("Data Distribution of DKO1 CGI Promoter UU") +
  labs(x = "log2(value)") + guides(color = FALSE) +
  scale_fill_discrete(labels = c("Rep1", "Rep2")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

# HCT116 non-CGI promoter MU
tmp_r1 <- melt(hct116_cpgpoor_MU$Replicate1)
tmp_r2 <- melt(hct116_cpgpoor_MU$Replicate2)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("Data Distribution of HCT116 Non-CGI promoter MU") +
  labs(x = "log2(value)") + guides(color = FALSE) +
  scale_fill_discrete(labels = c("Rep1", "Rep2")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

# HCT116 non-CGI promoter UU
tmp_r1 <- melt(hct116_cpgpoor_UU$Replicate1)
tmp_r2 <- melt(hct116_cpgpoor_UU$Replicate2)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("Data Distribution of HCT116 Non-CGI Promoter UU") +
  labs(x = "log2(value)") + guides(color = FALSE) +
  scale_fill_discrete(labels = c("Rep1", "Rep2")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

# DKO1 non-CGI promoter MU
tmp_r1 <- melt(dko1_cpgpoor_MU$Replicate1)
tmp_r2 <- melt(dko1_cpgpoor_MU$Replicate2)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("Data Distribution of DKO1 Non-CGI Promoter MU") +
  labs(x = "log2(value)") + guides(color = FALSE) +
  scale_fill_discrete(labels = c("Rep1", "Rep2")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

# DKO1 non-CGI promoter UU
tmp_r1 <- melt(dko1_cpgpoor_UU$Replicate1)
tmp_r2 <- melt(dko1_cpgpoor_UU$Replicate2)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  facet_wrap(~variable, scales = "free") +
  ggtitle("Data Distribution of DKO1 Non-CGI Promoter UU") +
  labs(x = "log2(value)") + guides(color = FALSE) +
  scale_fill_discrete(labels = c("Rep1", "Rep2")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

# H3K9me3 in CGI promoter MU (Rep1)
tmp_r1 <- melt(hct116_cpgrich_MU$Replicate1$H3K9me3)
tmp_r2 <- melt(dko1_cpgrich_MU$Replicate1$H3K9me3)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  ggtitle("Data Distribution of H3K9me3 in CGI Promoter MU (Rep1)") +
  labs(x = "log2(value)") + guides(color = "none") +
  scale_fill_discrete(labels = c("HCT116", "DKO1")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

# H3K9me3 in CGI promoter MU (Rep2)
tmp_r1 <- melt(hct116_cpgrich_MU$Replicate2$H3K9me3)
tmp_r2 <- melt(dko1_cpgrich_MU$Replicate2$H3K9me3)

ggplot(tmp_r1, aes(x = value, color = "#8DD3C7", fill = "#8DD3C7")) +
  geom_density(alpha = 0.7, size = 0.7) +
  geom_density(data = tmp_r2, alpha = 0.5, size = 0.7, aes(color = "#FFFFB3", fill = "#FFFFB3")) +
  ggtitle("Data Distribution of H3K9me3 in CGI Promoter MU (Rep2)") +
  labs(x = "log2(value)") + guides(color = "none") +
  scale_fill_discrete(labels = c("HCT116", "DKO1")) +
  theme_bw() + theme(legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     panel.grid = element_blank())

################################################################################
# Generate sactterplot matrix
################################################################################

i <- 1
j <- 1
dat <- summarized_matrix_list[[i]][[r[j]]]
ggpairs(dat, title = paste0("Scatterplot matrix of ", names(summarized_matrix_list)[i], "_", r[j]),
        upper = list(continuous = wrap("cor", method = "spearman")))
