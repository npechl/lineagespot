perc_of_lineages <- function(lineages_month){
  q <- lineages_month
  df <- as.data.frame(cbind(vec = unique(q), n = tabulate(match(q, unique(q)))))
  colnames(df) <- c('Lineage', 'Percentage')
  df$Percentage <- as.numeric(df$Percentage)
  df$Percentage <- df$Percentage/sum(df$Percentage)
  df <- df[order(df$Percentage,decreasing = T),]
  return(df)
}

read_xl_many_sheets <- function(path){
  sheetnames <- excel_sheets(path)
  lineages_pred <- lapply(excel_sheets(path), read_excel, path = path)
  for (i in 1:length(lineages_pred)){
    lineages_pred[[i]] <- as.data.frame(lineages_pred[[i]])
  }
  names(lineages_pred) <- sheetnames
  return(lineages_pred)
}

collapse_lineages_table <- function(lineages_pred){
  
  lineages <- unique(lineages_pred$Lineage)
  
  lineages_new <- NULL
  
  for (lin in lineages){
    lin_positions <- which(lineages_pred$Lineage == lin)
    lin_table <- lineages_pred[lin_positions, ]
    
    max_tree_ratio <- max(lin_table$Tree.Ratio)
    max_total_ratio <- max(lin_table$Total.Ratio)
    max_tree_av_ad <- max(lin_table$Tree.Avg.AD)
    max_total_av_ad<- max(lin_table$Total.Avg.AD)
    
    tree_av_av_ad <- mean(lin_table$Tree.Avg.AD)
    total_av_av_ad <- mean(lin_table$Total.Avg.AD)
    
    lineages_new <- rbind(lineages_new, c(lin,max_tree_ratio,max_total_ratio, max_tree_av_ad, max_total_av_ad, tree_av_av_ad, total_av_av_ad))
  }

  lineages_new <- as.data.frame(lineages_new)
  colnames(lineages_new) <- c("Lineage", "Max.Tree.Ratio","Max.Total.Ratio", "Max.Tree.Av.AD", "Max.Total.Av.AD", "Tree.Av.Av.AD", "Total.Av.Av.AD")
  
  lineages_new <- lineages_new[order(lineages_new$Max.Total.Ratio, decreasing = T),]
  
  row.names(lineages_new) <- c(1:dim(lineages_new)[1])
  
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
    
    pos <- pos + 1
  }
  
  return(lineages_pred)
}