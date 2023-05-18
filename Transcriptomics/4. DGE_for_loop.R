# Calculate average Log2FC of each DEG gene list associated with each metascape term
# Code written by Caroline He
# 2023

#STEP 1: SET UP ENVIORNMENT
library("dplyr")

setwd("~/multiome/CADASIL_transcriptomics_ONLY")

# 1. Specify DEG direction and celltype, load data
Celltype<- "Endo"
DGE <- 'UP'

## dataframe saved from celltype_ID_and_DEG script
DGE_df <- read_csv(paste(print(Celltype),"_CADASILvsControl_DESeq2.csv", sep = ""))

## term_gene is a dataframe with each column labeled as the metascape term ID, and each column's cells is populated by the genes associated with that term
term_genes <- read_csv(paste("metascape/", print(DGE), "_", print(Celltype), "_Metascape_genes.csv", sep = "")) 
## term_stats is a dataframe, with the metascape statistics for each term, it contains the following columns
## Term (metascape term ID), Descriptio (given description of term), metascape LogP (given p-value associated with term from metascape), metascape LogQ (given p-adj value associated with term from metascape), genes (genes from DEG list found to be associated with term)
term_stats <- read_csv(paste("metascape/", print(DGE), "_", print(Celltype), "_Metascape_stats.csv", sep = ""))

## 2. Create an empty dataframe to store final data
summary_cols = c('Term', 'Description', 'mean Log2FC', 'metascape LogP', 'metascape LogQ', 'genes')
summary_df = data.frame(matrix(nrow = 0, ncol = length(summary_cols))) 
colnames(summary_df) = summary_cols

## 3. Run a for loop to calculate the average Log2FC threshold associated with the DEGs found to be associated with a metascape term
for(i in 1:length(rownames(term_stats))){
  summary_df[i, 'Term'] = term_stats[i, 'Term']
  summary_df[i, 'Description'] = term_stats[i, 'Description']
  summary_df[i, 'metascape LogP'] = term_stats[i, 'LogP']
  summary_df[i, 'metascape LogQ'] = term_stats[i, 'Log(q-value)']
  summary_df[i, "genes"] = term_stats[i,'Symbols']
  gene_list = pull(term_genes, i)
  gene_list_DGEs = DGE_df[DGE_df$"...1" %in% gene_list,]
  gene_list_DGEs = pull(gene_list_DGEs, 'avg_log2FC')
  term_avg_DGE = mean(gene_list_DGEs)
  summary_df[i, 'mean Log2FC'] = term_avg_DGE
  
}

write.csv(summary_df, paste(print(DGE), "_", print(Celltype), "_mean_DGE_per_term.csv", sep = ""))

