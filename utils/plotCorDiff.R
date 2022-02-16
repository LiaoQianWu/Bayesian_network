library(reshape2)
library(ggplot2)
library(ggrepel)

path <- "/path/to/summarizedMatrix/"

hct116_cpgrich_MU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgrich_MU_allEpigenomes.RDS"))
hct116_cpgrich_UU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgrich_UU_allEpigenomes.RDS"))
hct116_cpgpoor_MU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgpoor_MU_allEpigenomes.RDS"))
hct116_cpgpoor_UU <- readRDS(paste0(path, "All_Consolidated_HCT116/hct116_cpgpoor_UU_allEpigenomes.RDS"))
dko1_cpgrich_MU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgrich_MU_allEpigenomes.RDS"))
dko1_cpgrich_UU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgrich_UU_allEpigenomes.RDS"))
dko1_cpgpoor_MU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgpoor_MU_allEpigenomes.RDS"))
dko1_cpgpoor_UU <- readRDS(paste0(path, "All_Consolidated_DKO1/dko1_cpgpoor_UU_allEpigenomes.RDS"))

################################################################################
# Compute correlation difference between MU and UU
################################################################################

corDiff <- function(MU, UU){
  avg_MU <- (MU$Replicate1 + MU$Replicate2) / 2
  avg_UU <- (UU$Replicate1 + UU$Replicate2) / 2
  
  if(sum(avg_MU$Methylation) == 0 | sum(avg_UU$Methylation) == 0){
    avg_MU <- avg_MU[, -which(colnames(avg_MU) == "Methylation")]
    avg_UU <- avg_UU[, -which(colnames(avg_UU) == "Methylation")]
  }
  
  cor_MU <- cor(avg_MU)
  cor_UU <- cor(avg_UU)
  
  cor_melt_MU <- melt(cor_MU)
  cor_melt_UU <- melt(cor_UU)
  
  id_MU <- paste(cor_melt_MU$Var1, "&", cor_melt_MU$Var2)
  id_UU <- paste(cor_melt_UU$Var1, "&", cor_melt_UU$Var2)
  
  if(identical(id_MU, id_UU)){
    cor_mat <- cbind(cor_melt_MU, cor_melt_UU$value)
    rownames(cor_mat) <- id_MU
    colnames(cor_mat) = c("Feature1", "Feature2", "cor_MU", "cor_UU")
    cor_mat$diff <- cor_mat$cor_MU - cor_mat$cor_UU
    cor_mat <- cor_mat[order(cor_mat$diff),]
    itself <- cor_mat[which(cor_mat$diff == 0),]
    tmp <- cor_mat[which(cor_mat$diff != 0),]
    tmp <- tmp[seq(1, nrow(tmp), 2),]
    cor_mat <- rbind(tmp, itself)
    cor_mat <- cor_mat[order(cor_mat$diff),]
    cor_mat$avg_cor <- (abs(cor_mat$cor_MU) + abs(cor_mat$cor_UU)) / 2
    cor_mat$id <- 1:nrow(cor_mat)
    cor_mat <- cor_mat[, c("id", "diff", "avg_cor")]
  }
  return(cor_mat)
}

################################################################################
# Plot correlation difference between HCT116 and DKO1
################################################################################

# HCT116 and DKO1 CGI promoter MU
diff_cpgrich_MU <- corDiff(hct116_cpgrich_MU, dko1_cpgrich_MU)
POI <- diff_cpgrich_MU[diff_cpgrich_MU$diff > 0.2 | diff_cpgrich_MU$diff < -0.2,]

ggplot(diff_cpgrich_MU, aes(x = id, y = diff, size = avg_cor)) +
  geom_point() +
  geom_point(data = POI, aes(y = diff), color = "red", show.legend = F) +
  ggtitle("Difference of Correlation Between HCT1116 and DKO1 CGI Promoter MU") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()) +
  scale_fill_manual(values = "black") +
  labs(x = "index", y = "difference") +
  ylim(-0.4, 0.3) +
  geom_text_repel(data = POI, aes(label = rownames(POI)), size = 3.3)

# HCT116 and DKO1 CGI non-promoter MU
diff_cpgpoor_MU <- corDiff(hct116_cpgpoor_MU, dko1_cpgpoor_MU)
POI <- diff_cpgpoor_MU[diff_cpgpoor_MU$diff > 0.1 | diff_cpgpoor_MU$diff < -0.1,]

ggplot(diff_cpgpoor_MU, aes(x = id, y = diff, size = avg_cor)) +
  geom_point() +
  geom_point(data = POI, aes(y = diff), color = "red", show.legend = F) +
  ggtitle("Difference of Correlation Between HCT1116 and DKO1 Non-CGI Promoter MU") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()) +
  scale_fill_manual(values = "black") +
  labs(x = "index", y = "difference") +
  ylim(-0.5, 0.2) +
  geom_text_repel(data = POI, aes(label = rownames(POI)), size = 3.3)
