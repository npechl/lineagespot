# Libraries needed ----------------------------------

library(data.table)
library(tidyverse)
library(saveImageHigh)
library(ComplexHeatmap)
library(circlize)
library(gdata)
library(stringdist)
library(ggrepel)
library(hrbrthemes)
library(viridis)
library(xlsx)

base::rm(list = ls())

source("plot_functions.R")

# Inputs ----------------------------------
# mean_all --> mean of all the existing mutations
# simple_mean --> mean of  all the mutations
# min_existing_unique --> min of the unique existing (greater than 0) mutation
# mean_existing_unique --> mean of the existing unique mutations
# simple_mean_unique --> mean of the unique mutations

method <- "mean_existing_unique"

input_data <- read.csv("new files/sars-sewage.csv")

# Read meta-data file ------------

meta_file <- "new files/meta-data.xlsx"
meta_data <- read.xlsx(meta_file, sheetIndex = 1)

timepoints <- data.table(samples = meta_data$Sample_ID,
                         heat_labels = meta_data$Time.length_.sampling.)

rm(meta_data, meta_file)

# Read amino acids ------------------

amino_acids <- "Amino_acid_Abbreviations.txt"
aa_abbreviations <- read.table(amino_acids, header = T, sep = ",")

# Set mutations source and filepaths ----------------

mutations_source <-  "outbreak_info" # "pangolin"   #or "veo" 

plots_folder <- "new sars/new samples/outbreak_info"

files <- list.files("outbreak_info/")

# Create heatmaps ------------------------

input_data <- input_data[which(input_data$ID != ""), ]
col_fun = colorRamp2(c( 0,1), c("#FFF6FF", "red"))

all_data <- list()

if (mutations_source == "outbreak_info"){
  for (specific_fp in files){
    
    specific_fp <- paste0("outbreak_info/", specific_fp)
    aa_abbreviations <- read.table(amino_acids, header = T, sep = ",")
    
    specific_lineage <- read.table(specific_fp, header = T, sep = "\t")
    
    deletions <- specific_lineage[which(specific_lineage$type == "deletion"), ]
    others <- specific_lineage[which(specific_lineage$type == "substitution"), ]
    
    others$ref_aa <-  paste0(others$ref_aa, " ")
    
    for (i in c(1:nrow(aa_abbreviations))){
      
      others$ref_aa <- str_replace_all(others$ref_aa, 
                                       paste0(aa_abbreviations$One_Letter[i], " "), 
                                       aa_abbreviations$Three_Letter[i])
    }
    
    others$alt_aa <-  paste0(others$alt_aa, " ")
    
    for (i in c(1:nrow(aa_abbreviations))){
      
      others$alt_aa <- str_replace_all(others$alt_aa,
                                       paste0(aa_abbreviations$One_Letter[i]," "), 
                                       aa_abbreviations$Three_Letter[i])
    }
    
    
    
    list_aa <- check_amino(input_data,
                           others,
                           timepoints)
    
    if (nrow(deletions) > 0){
      
      list_del <- check_deletions(input_data,
                                  deletions,
                                  timepoints,
                                  aa_abbreviations)
      
      heat_matrix <- rbind(list_del$data_specific, list_aa$data_specific)
      
      del_annot <- list_del$annotations
      annotations <- rbind(del_annot[, -c("match")], list_aa$annotations)
      row.names(annotations) <- c(del_annot$match, row.names(list_aa$annotations))
      
      
    } else {
      
      heat_matrix <- list_aa$data_specific
      annotations <- list_aa$annotations
    }
    
    
    # 
    # row.names(annotations) <- annotations$mutation
    # row.names(heat_matrix) <- annotations$mutation
    
    #heat_matrix <- heat_matrix[specific_lineage$mutation, ]
    #annotations <- annotations[c(row.names(heat_matrix)), ]
    #annotations <- annotations[match(specific_lineage$mutations$V1, row.names(annotations)),]
    annotations$names <- row.names(annotations)
    annotations$names <- reorder.factor(annotations$names, new.order=row.names(heat_matrix))
    new_annotations <- annotations %>%
      arrange(names)
    
    new_annotations[is.na(new_annotations)] <- " "
    new_annotations[new_annotations == ""] <- " "
    
    
    specific_fp <- str_remove(specific_fp, ".tsv")
    specific_fp <- str_remove(specific_fp, "outbreak_info/outbreakinfo_mutation_report_data_")
    
    
    all_data[[specific_fp]] <- cbind(heat_matrix, new_annotations$mutation)
  } 
  
} else if (mutations_source == "pangolin"){
  
  for (specific_fp in files){
    
    specific_fp <- paste0("pangolin data/", specific_fp)
    aa_abbreviations <- read.table(amino_acids, header = T, sep = ",")
    
    specific_lineage <- read.table(specific_fp, header = F, sep = "\t")
    
    specific_lineage <- specific_lineage %>%
      separate(V1, c("type", "gene", "HGVS.p"), sep =":", remove = FALSE)
    
    colnames(specific_lineage) <- c("label", "type", "gene", "HGVS.p")
    
    others <- specific_lineage[which(specific_lineage$type == "aa"), ]
    
    for (i in c(1:nrow(aa_abbreviations))){
      
      others$HGVS.p <- str_replace_all(others$HGVS.p, 
                                       aa_abbreviations$One_Letter[i], 
                                       paste0(aa_abbreviations$One_Letter[i], " "))
    }
    
    for (i in c(1:nrow(aa_abbreviations))){
      
      others$HGVS.p <- str_replace_all(others$HGVS.p, 
                                       paste0(aa_abbreviations$One_Letter[i], " "), 
                                       paste0(aa_abbreviations$Three_Letter[i], " "))
    }
    
    others$HGVS.p <-  str_remove_all(others$HGVS.p, " ")
    
    list_aa <- pangolin_amino(input_data,
                              others,
                              timepoints)
    
    deletions <- specific_lineage[which(specific_lineage$type == "del"), ]
    
    if (nrow(deletions) > 0){
      
      colnames(deletions) <- c("label", "type", "position", "length")
      
      list_del <- pangolin_deletions(input_data,
                                     deletions,
                                     timepoints,
                                     aa_abbreviations)
      
      heat_matrix <- rbind(list_del$data_specific, list_aa$data_specific)
      
      annotations <- rbind(list_del$annotations, list_aa$annotations)
      
      row.names(annotations) <- annotations$label
      
    } else {
      
      heat_matrix <- list_aa$data_specific
      annotations <- list_aa$annotations
    }
    
    
    # 
    # row.names(annotations) <- annotations$mutation
    # row.names(heat_matrix) <- annotations$mutation
    
    #heat_matrix <- heat_matrix[specific_lineage$mutation, ]
    #annotations <- annotations[c(row.names(heat_matrix)), ]
    
    #annotations <- annotations[match(specific_lineage$mutations$V1, row.names(annotations)),]
    #annotations$names <- row.names(annotations)
    
    row.names(heat_matrix) <- row.names(annotations)
    heat_matrix <- heat_matrix[specific_lineage$label, ]
    
    annotations$label <- reorder.factor(annotations$label, new.order=row.names(heat_matrix))
    new_annotations <- annotations %>%
      arrange(label)
    
    new_annotations[is.na(new_annotations)] <- " "
    new_annotations[new_annotations == ""] <- " "
    
    #row.names(heat_matrix) <- new_annotations$mutation
    specific_fp <- str_remove(specific_fp, ".txt")
    specific_fp <- unlist(str_split(specific_fp, "/"))
    specific_fp <- specific_fp[length(specific_fp)]
    
    new_annotations$mutation <- new_annotations$label
    
    all_data[[specific_fp]] <- cbind(heat_matrix, row.names(new_annotations))
  }
  
} else if (mutations_source == "veo"){
  
  for (specific_fp in files){
    
    specific_fp <- paste0("veo table2 data/", specific_fp)
  
    specific_lineage <- read.table(specific_fp, header = F, sep = "\t")
    
    specific_lineage <- specific_lineage %>%
      separate(V1, c("type", "gene", "HGVS.p"), sep =":", remove = FALSE)
    
    colnames(specific_lineage) <- c("label", "type", "gene", "HGVS.p")
    
    for (i in c(1:nrow(aa_abbreviations))){
      
      specific_lineage$HGVS.p <- str_replace_all(specific_lineage$HGVS.p, 
                                                 aa_abbreviations$One_Letter[i], 
                                                 paste0(aa_abbreviations$One_Letter[i], " "))
    }
    
    for (i in c(1:nrow(aa_abbreviations))){
      
      specific_lineage$HGVS.p <- str_replace_all(specific_lineage$HGVS.p, 
                                                 paste0(aa_abbreviations$One_Letter[i], " "), 
                                                 paste0(aa_abbreviations$Three_Letter[i], " "))
    }
    
    specific_lineage$HGVS.p <-  str_remove_all(specific_lineage$HGVS.p, " ")
    
    others <- specific_lineage[which(specific_lineage$type == "aa"), ]
    
    list_aa <- pangolin_amino(input_data,
                              others,
                              timepoints)
    
    deletions <- specific_lineage[which(specific_lineage$type == "del"), ]
    
    if (nrow(deletions) > 0){
      
      colnames(deletions) <- c("label", "type", "gene", "HGVS.p")
      
      list_del <- veo_deletions(input_data,
                                deletions,
                                timepoints,
                                aa_abbreviations)
      
      heat_matrix <- rbind(list_del$data_specific, list_aa$data_specific)
      
      annotations <- rbind(list_del$annotations, list_aa$annotations)
      
      row.names(annotations) <- annotations$label
      
    } else {
      
      heat_matrix <- list_aa$data_specific
      annotations <- list_aa$annotations
    }
    
    
    # 
    # row.names(annotations) <- annotations$mutation
    # row.names(heat_matrix) <- annotations$mutation
    
    #heat_matrix <- heat_matrix[specific_lineage$mutation, ]
    #annotations <- annotations[c(row.names(heat_matrix)), ]
    
    #annotations <- annotations[match(specific_lineage$mutations$V1, row.names(annotations)),]
    #annotations$names <- row.names(annotations)
    
    row.names(heat_matrix) <- row.names(annotations)
    heat_matrix <- heat_matrix[specific_lineage$label, ]
    
    annotations$label <- reorder.factor(annotations$label, new.order=row.names(heat_matrix))
    new_annotations <- annotations %>%
      arrange(label)
    
    new_annotations[is.na(new_annotations)] <- " "
    new_annotations[new_annotations == ""] <- " "
    
    #row.names(heat_matrix) <- new_annotations$mutation
    specific_fp <- str_remove(specific_fp, ".txt")
    specific_fp <- unlist(str_split(specific_fp, "/"))
    specific_fp <- specific_fp[length(specific_fp)]
    
    new_annotations$mutation <- new_annotations$label
    
    all_data[[specific_fp]] <- cbind(heat_matrix, row.names(new_annotations))
    
  } 
}


unique_mutations <- list()

if ((method != "mean_all") & (method != "simple_mean")) { 
  
  for (j in c(1:length(all_data))){
    
    iterations <- setdiff(c(1:length(all_data)), j)
    
    temp <- all_data[[j]]
    temp <- temp[, nrow(timepoints) +1]
    
    for (i in iterations){
      
      sth_else <- all_data[[i]]
      sth_else <- sth_else[, nrow(timepoints) +1]
      temp <- setdiff(temp, sth_else)
      
    }
    
    unique_mutations[[j]] <- temp
  }
  
} else {
  
  for (j in c(1:length(all_data))) {
    
    temp <- all_data[[j]]
    unique_mutations[[j]] <- temp[, nrow(timepoints) +1]
    
  }
}


table_means <- data.table(lineage = character(),
                          December_02_14 = numeric(),
                          February_05_11 = numeric(),
                          February_12_18 = numeric(),
                          February_19_25 = numeric(),
                          February_26_March_04 = numeric(),
                          March_05_11 = numeric(),
                          March_12_18 = numeric(),
                          March_19_25 = numeric(),
                          March_26_April_01 = numeric(),
                          April_02_08 = numeric(),
                          April_09_15 = numeric())


  
for (i in c(1:length(all_data))){
    
    temp <- all_data[[i]]
    temp <- temp[which(temp[, nrow(timepoints) +1] %in% unique_mutations[[i]]), ]
    
    if (str_detect(method, "simple", negate = TRUE)){
      temp[temp == "0"] <- NA
    }
    
    if (is.vector(temp)) {
      temp[is.na(temp)] <- 0
      one_row <- data.table(lineage = names(all_data)[i],
                            December_02_14 = temp[1],
                            February_05_11 =  temp[2],
                            February_12_18 =  temp[3],
                            February_19_25 = temp[4],
                            February_26_March_04 = temp[5],
                            March_05_11 =  temp[6],
                            March_12_18 = temp[7],
                            March_19_25 = temp[8],
                            March_26_April_01 = temp[9],
                            April_02_08 = temp[10],
                            April_09_15 = temp[11])
      
    } else if (nrow(temp) == 0) {
      temp[is.na(temp)] <- 0
      one_row <- data.table(lineage = names(all_data)[i],
                            December_02_14 = temp[1],
                            February_05_11 =  temp[2],
                            February_12_18 =  temp[3],
                            February_19_25 = temp[4],
                            February_26_March_04 = temp[5],
                            March_05_11 =  temp[6],
                            March_12_18 = temp[7],
                            March_19_25 = temp[8],
                            March_26_April_01 = temp[9],
                            April_02_08 = temp[10],
                            April_09_15 = temp[11])
      
    } else if (is.matrix(temp)){
      
      temp <- temp[, 1:11]
      
      if (is.vector(temp)) {
        temp[is.na(temp)] <- 0
        
        one_row <- data.table(lineage = names(all_data)[i],
                              December_02_14 = vector_[1],
                              February_05_11 =  vector_[2],
                              February_12_18 =  vector_[3],
                              February_19_25 = vector_[4],
                              February_26_March_04 = vector_[5],
                              March_05_11 =  vector_[6],
                              March_12_18 = vector_[7],
                              March_19_25 = vector_[8],
                              March_26_April_01 = vector_[9],
                              April_02_08 =vector_[10],
                              April_09_15 = vector_[11]) 
      } else {
        
        temp <- apply(temp, 2, as.numeric)
        
        if (str_detect(method, "mean")) {
          
          vector_ <- colMeans(temp, na.rm = T)
          vector_[is.na(vector_)] <- 0
          
        } else if (str_detect(method, "min")) {
          
          vector_ <- numeric(11)
          for (x in 1:ncol(temp)) {
            vector_[x] <- min(temp[, x], na.rm = TRUE)
          }
          
          vector_[vector_ == "Inf"] <- 0
        }
        
        
        
        one_row <- data.table(lineage = names(all_data)[i],
                              December_02_14 = vector_[1],
                              February_05_11 =  vector_[2],
                              February_12_18 =  vector_[3],
                              February_19_25 = vector_[4],
                              February_26_March_04 = vector_[5],
                              March_05_11 =  vector_[6],
                              March_12_18 = vector_[7],
                              March_19_25 = vector_[8],
                              March_26_April_01 = vector_[9],
                              April_02_08 =vector_[10],
                              April_09_15 = vector_[11]) 
      }
      
    }
    
    #one_row[1, 2:12] <- one_row[1, 2:12] * 100
    table_means <- rbind(table_means, one_row) 
  }


colnames(table_means)[2:ncol(table_means)] <- timepoints$heat_labels

write.csv(table_means, paste0(plots_folder,"/", method, "_vector_values.csv"), row.names = F)

