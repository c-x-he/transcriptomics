#Statistical tests
library(tidyverse)
library(rstatix)
library(ggpubr)
library(effsize)  
library(ggplot2)

setwd("~/TREM2/hdWGCNA/Oli-subset")
seurat_obj <- readRDS("~/TREM2/hdWGCNA/Oli-subset/Oli_hdWGCNA-seurat_obj-basic-end.rds")

modules <- c("Oli-M1", "Oli-M2", "Oli-M3", "Oli-M4", "Oli-M5", "Oli-M6", "Oli-M7", "Oli-M8")
modules <- data.frame(modules)

#########Kruskal_Wallis test#########
Kruskal_Wallis<- setNames(data.frame(matrix(ncol = 3, nrow = nrow(modules))), c("Module", "Kruskal-Wallis chi-squared", "p-value"))
for (i in 1:nrow(modules)){
  hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[[paste(modules[i,1])]])
  rownames(hME)=rownames(seurat_obj@meta.data)
  
  #Kruskal-Wallis test
  Kruskal_Wallis_data<-kruskal.test(hME ~ diagnosis, data = hME)
  
  Kruskal_Wallis[i,1]<-modules[i,1]
  Kruskal_Wallis[i,2]<-Kruskal_Wallis_data$statistic
  Kruskal_Wallis[i,3]<-Kruskal_Wallis_data$p.value
}

write_csv(Kruskal_Wallis, "Kruskal_Wallis_data.csv")


#########Wilcoxon Rank Sum test#########
#TREM2 vs Control  
Wilcoxon_rank_sum<- setNames(data.frame(matrix(ncol = 5, nrow = nrow(modules))), c("Module", "Wilcoxon p-value (BH adjusted)", "Delaney's A effect size", "Delaney's A magnitude", "Wilcoxon effect size"))

for (i in 1:nrow(modules)){
  hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[[paste(modules[i,1])]])
  rownames(hME)=rownames(seurat_obj@meta.data)
  
  # Show a sample of the data by group
  set.seed(223)
  
  #signifiance test
  stat.test <- hME %>% 
    rstatix::wilcox_test(hME ~ diagnosis) %>%
    add_significance()
  stat.test
  W_effect_size<-hME %>% wilcox_effsize(hME ~ diagnosis)
  
  Wilcoxon_rank_sum[i,1] <- modules[i,1]
  Wilcoxon_rank_sum[i,2] <- stat.test[1,7]
  Wilcoxon_rank_sum[i,5] <- W_effect_size[1,4]
         
  #get AUC effect size
  TREM2<-hME$hME[hME$diagnosis=='TREM2 Homozy']
  Control<-hME$hME[hME$diagnosis=='Control']
         
  Delaneys_A<-VD.A(TREM2, Control)
  Wilcoxon_rank_sum[i,3] <-Delaneys_A[["estimate"]]
  Wilcoxon_rank_sum[i,4] <-Delaneys_A[["magnitude"]]   
}

write_csv(Wilcoxon_rank_sum, "Wilcoxon_rank_sum_data.csv")


########IN PROGRESS: DO NOT USE IM LOOKING FOR A BUG########################
#########Graphs for Wilcoxon#########
for (i in 1:nrow(modules)){
  
  hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[[paste(modules[i,1])]])
  rownames(hME)=rownames(seurat_obj@meta.data)
  
  #make boxplot
  bxp <-ggboxplot(
    hME, x = "diagnosis", y = "hME", 
    ylab = "hME", xlab = "Diganosis", add = "jitter"
  )
  
  #add to graph
  stat.test <- stat.test %>% add_xy_position(x = "diagnosis")
  bxp <-bxp + 
    stat_pvalue_manual(stat.test, tip.length = 0) +
    labs(title = paste(modules[i,1]))
  assign(paste('bxp_', print(modules[i,1]),sep = ""), bxp)
}


`bxp_Oli-M1`|`bxp_Oli-M2`|`bxp_Oli-M3`|`bxp_Oli-M4`
`bxp_Oli-M5`|`bxp_Oli-M6`|`bxp_Oli-M7`|`bxp_Oli-M8`



