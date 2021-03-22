rm(list = ls())

library(openxlsx)

# PARAMETERS ----------------------------

# Specify thresholds for tree ration and total ratio
# total_ratio_threshold <- 0.75
# tree_ratio_threshold <- 0 

# INPUT FILES ---------------------------

# File to collapse
analysis_output_file <- 'vcf-files/Sewage_CoV19_L3_S1_L001_freebayes.ann-lineagespot.tsv'
# analysis_output_folder <- paste(dirname(getwd()), 'analysis-output', sep = '/')
# analysis_output_fullpath <- paste(analysis_output_folder, analysis_output_file, sep = '/')

# ----------------- Functions --------------------------------------------------
collapse_lineages_table <- function(lineages_pred){
  
  lineages <- unique(lineages_pred$Lineage)
  
  lineages_new <- list()
  
  for (lin in lineages){
    
    # lin_positions <- which(lineages_pred$Lineage == lin)
    lin_table <- lineages_pred[which(lineages_pred$Lineage == lin), ]
    
    # pos_of_max <- which(lin_table$Total.Ratio == max(lin_table$Total.Ratio))
    
    mean_tree_ratio <- mean(lin_table$Tree.Ratio)
    mean_total_ratio <- mean(lin_table$Total.Ratio)
    mean_total_ratio_var <- mean(lin_table$Total.Ratio.Var)
    mean_tree_av_ad <- mean(lin_table$Tree.Avg.AD)
    mean_total_av_ad<- mean(lin_table$Total.Avg.AD)
    mean_sig = mean(lin_table$Sig)
    
    lineages_new[[lin]] <- data.table::data.table("Lineage" = lin, 
                                                  "Mean.Tree.Ratio" = mean_tree_ratio,
                                                  "Mean.Total.Ratio" = mean_total_ratio, 
                                                  "Mean.Total.Ratio.Var" = mean_total_ratio_var,
                                                  "Mean.Tree.Av.AD" = mean_tree_av_ad, 
                                                  "Mean.Total.Av.AD" = mean_total_av_ad,
                                                  "Mean.Sig" = mean_sig)
  }
  
  # lineages_new <- as.data.frame(lineages_new)
  # colnames(lineages_new) <- c()
  
  lineages_new = data.table::rbindlist(lineages_new)
  
  # lineages_new <- lineages_new[order(lineages_new$Max.Total.Ratio, decreasing = TRUE),]
  
  # row.names(lineages_new) <- 1:nrow(lineages_new)
  
  return(lineages_new)
  
}

add_columns_to_initial_matrix <- function(lineages_pred, lineages_collapsed){
  
  lineages_pred$Max.Tree.Ratio <- 0
  lineages_pred$Max.Total.Ratio <- 0
  lineages_pred$Max.Tree.Av.AD <- 0
  lineages_pred$Max.Total.Av.AD <- 0
  lineages_pred$Tree.Av.Av.AD <- 0
  lineages_pred$Total.Av.Av.AD <- 0
  
  lins <- lineages_collapsed$Lineage
  pos <- 1
  for (lin in lins){
    lin_positions <- which(lineages_pred$Lineage == lin)
    
    lineages_pred[lin_positions,]$Max.Tree.Ratio <- lineages_collapsed[pos,]$Max.Tree.Ratio
    lineages_pred[lin_positions,]$Max.Total.Ratio <- lineages_collapsed[pos,]$Max.Total.Ratio
    lineages_pred[lin_positions,]$Max.Tree.Av.AD <- lineages_collapsed[pos,]$Max.Tree.Av.AD
    lineages_pred[lin_positions,]$Max.Total.Av.AD <- lineages_collapsed[pos,]$Max.Total.Av.AD
    lineages_pred[lin_positions,]$Tree.Av.Av.AD <- lineages_collapsed[pos,]$Tree.Av.Av.AD
    lineages_pred[lin_positions,]$Total.Av.Av.AD <- lineages_collapsed[pos,]$Total.Av.Av.AD
    
    # Change total avg ad
    pos <- pos + 1
  }
  
  return(lineages_pred)
}

# Importing results from our analysis -----------------------


# Import
lineages_pred <- read.csv(analysis_output_file, sep = '\t')

# Filter lineages based on thresholds
lineages_pred <- lineages_pred[which(lineages_pred$Total.Overlap.Var != 0), ]
# lineages_pred <- lineages_pred[which(lineages_pred$Total.Ratio >= total_ratio_threshold), ]
# lineages_pred <- lineages_pred[which(lineages_pred$Tree.Ratio >= tree_ratio_threshold), ]

# Collapse
lineages_collapsed <- collapse_lineages_table(lineages_pred)

# Add new columns to initial table
lineages_pred <- add_columns_to_initial_matrix(lineages_pred, lineages_collapsed )

# Save tables
utils::write.table(lineages_pred, 
                   file = stringr::str_replace(analysis_output_file, ".tsv", "_modified.tsv"), 
                   row.names = FALSE, 
                   quote = FALSE, 
                   sep = "\t")

utils::write.table(lineages_collapsed, 
                   file = stringr::str_replace(analysis_output_file, ".tsv", "-collapsed.tsv"), 
                   row.names = FALSE, 
                   quote = FALSE, 
                   sep = "\t")
