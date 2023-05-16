#Example for loop for Eliza, possible downstream analysis after metascape + DGE, or hdWGCNA
library("dplyr")

setwd("~/new_frontal_official/hdWGCNA/Ast-subset")

module <- 'Ast-M2'
kME_df <- read_csv("Ast_modules.subset-first.csv")
term_genes <- read_csv(paste("Metascape/", print(module), "_Metascape_genes.csv", sep = ""))
term_stats <- read_csv(paste("Metascape/", print(module), "_Metascape_stats.csv", sep = ""))

#create an empty dataframe to store final data
summary_cols = c('Term', 'Description', 'mean kME', 'metascape LogP', 'metascape LogQ')
summary_df = data.frame(matrix(nrow = 0, ncol = length(summary_cols))) 
colnames(summary_df) = summary_cols

for(i in 1:length(rownames(term_stats))){
  summary_df[i, 'Term'] = term_stats[i, 'Term']
  summary_df[i, 'Description'] = term_stats[i, 'Description']
  summary_df[i, 'metascape LogP'] = term_stats[i, 'LogP']
  summary_df[i, 'metascape LogQ'] = term_stats[i, 'Log(q-value)']
  gene_list = pull(term_genes, i)
  gene_list_kMEs = kME_df[kME_df$gene_name %in% gene_list,]
  module_colname<-paste('kME_', print(module), sep="")
  gene_list_kMEs = pull(gene_list_kMEs, module_colname)
  term_avg_kME = mean(gene_list_kMEs)
  summary_df[i, 'mean kME'] = term_avg_kME
}

write.csv(summary_df, paste(module, "_mean_kME_per_term.csv", sep = ""))

#graph after this


