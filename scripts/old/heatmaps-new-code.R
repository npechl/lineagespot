# Libraries needed ----------------------------------

library(data.table)
library(tidyverse)
library(saveImageHigh)
library(ComplexHeatmap)
library(circlize)
library(gdata)
library(stringdist)

base::rm(list = ls())

source("plot_functions.R")

# Inputs ----------------------------------

input_data <- read.csv("sars-sewage.csv")

timepoints <- data.table(samples = c("CoV19_L1_S1", "Sewage-L2_S10_L001", 
                                     "Sewage_CoV19_L3_S1_L001", "L4_S1_L001",             
                                     "L5_S2_L001", "L6_S3_L001", "L7_S4_L001"),
                         heat_labels = c("02-14/12/2020","05-11/02/2021", 
                                         "12-18/02/2021", "19-25/02/2021",
                                         "26/02-4/03/2021", "5-11/03/2021",
                                         "12-18/03/2021"))

specific_fp <- "outbreak_info/outbreakinfo_mutation_report_data_b.1.351.tsv"
amino_acids <- "outbreak_info/Amino_acid_Abbreviations.txt"

plots_folder <- "new sars/outbreak"

threshold <- 0.4

# Create heatmaps ------------------------

input_data <- input_data[which(input_data$ID != ""), ]
col_fun = colorRamp2(c( 0,1), c("#FFF6FF", "red"))

# Common variants --------------------------

#common variants
common_variants <- dplyr::count(input_data, ID)
common_variants <- common_variants[which(common_variants$n == nrow(timepoints)), ]
common_variants <- input_data[which(input_data$ID %in% common_variants$ID ), ]

common_variants$POS <- as.numeric(common_variants$POS)

common_variants <- common_variants[order(common_variants$POS, decreasing = F), ]

if (nrow(common_variants) > 0){
  annot_data <- unique(common_variants[, c("POS" ,"ID", "REF", "ALT", "Gene_Name", "HGVS.c", "HGVS.p")])
  
  annot <- paste0(annot_data$REF, " -> ", annot_data$ALT)
  
  genes <- annot_data$Gene_Name
  
  HGVS.c <- annot_data$HGVS.c
  
  HGVS.p <- annot_data$HGVS.p
  
  extra = rowAnnotation(variants = anno_text(annot, location = 0.5, just = "center",
                                             gp = gpar(fontsize = 11, fill = "#DCE5F3", col = "black", border = "white"),
                                             width = max_text_width(annot)*1.1),
                        Gene_name = anno_text(genes,
                                              location = 0.5, just = "center",
                                              gp = gpar(fontsize = 11, fill = "#F5CDF7", col = "black", border = "white"),
                                              width = max_text_width(genes)),
                        HGVS.c = anno_text(HGVS.c,
                                           location = 0.5, just = "center",
                                           gp = gpar(fontsize = 11, fill = "#D3F7CD", col = "black", border = "white"),
                                           width = max_text_width(HGVS.c)),
                        HGVS.p = anno_text(HGVS.p,
                                           location = 0.5, just = "center",
                                           gp = gpar(fontsize = 11, fill = "#F8E5A6", col = "black", border = "white"),
                                           width = max_text_width(HGVS.p))
  ) 
  
  data_common <- matrix(0, nrow = nrow(annot_data), ncol = nrow(timepoints))
  row.names(data_common) <- annot_data$POS
  colnames(data_common) <- timepoints$samples
  
  for (k in timepoints$samples){
    
    data_common[, k] <- common_variants[which(common_variants$sample == k), ]$AF
    
  }
  
  colnames(data_common) <- timepoints$heat_labels
  
  graph1 <-Heatmap(data_common, name = "frequency", 
                   col = col_fun,
                   heatmap_width = unit(18 +nrow(timepoints), "cm"), heatmap_height = unit(3+nrow(annot_data), "cm"),
                   border = T,
                   #border_gp = gpar(col = "grey"),
                   rect_gp = gpar(col = "white", lwd = 1),
                   cluster_rows = FALSE,cluster_columns  = FALSE,
                   column_names_gp = gpar(fontsize = 13),
                   row_names_gp = gpar(fontsize = 14),
                   column_names_rot = 50,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.2f", data_common[i,j]), x, y, gp = gpar(fontsize = 12))
                   },
                   heatmap_legend_param = list(
                     labels_gp = gpar(fontsize = 10),
                     title_gp =gpar(fontsize = 13) ,
                     legend_height = unit(2.5, "cm"),
                     legend_width = unit(4, "cm")
                   ),
                   right_annotation = extra,
                   row_split = factor(genes, levels = unique(genes)), #factor(rep(LETTERS[1:3], 6), levels = LETTERS[3:1])
                   row_order = order(as.numeric(row.names(data_common)), decreasing = F),row_title_rot = 0,
                   cluster_row_slices = F,
  )
  
  
  save_image(print(graph1),paste0(plots_folder, "/common variants.png"), 
             width = 10, height = 10, res = 100)  
}


# All variants --------------------------

#all variants

all_variants <- input_data[which(input_data$AF > threshold), ]

all_variants$POS <- as.numeric(all_variants$POS)

all_variants <-all_variants[order(all_variants$POS, decreasing = F), ]

annot_data <- unique(all_variants[, c("POS" ,"ID", "REF", "ALT", "Gene_Name", "HGVS.c", "HGVS.p")])

annot <- paste0(annot_data$REF, " -> ", annot_data$ALT)

who.huge <- which(nchar(annot) > 20)
annot[who.huge] <- paste0(str_sub(annot_data$REF[who.huge],1,10),"...", " -> ",
                          str_sub(annot_data$ALT[who.huge],1,10),"...")

genes <- annot_data$Gene_Name

HGVS.c <- annot_data$HGVS.c

who.huge <- which(nchar(HGVS.c) > 25)
HGVS.c[who.huge] <- paste0(str_sub(HGVS.c[who.huge],1,25),"...")

HGVS.p <- annot_data$HGVS.p

who.huge <- which(nchar(HGVS.p) > 25)
HGVS.p[who.huge] <- paste0(str_sub(HGVS.p[who.huge],1,25),"...")

extra = rowAnnotation(variants = anno_text(annot, location = 0.5, just = "center",
                                           gp = gpar(fontsize = 8, fill = "#DCE5F3", col = "black", border = "white"),
                                           width = unit(5, "cm")), # max_text_width(annot)),
                      Gene_name = anno_text(genes,
                                            location = 0.5, just = "center",
                                            gp = gpar(fontsize = 8, fill = "#F5CDF7", col = "black", border = "white"),
                                            width = unit(2, "cm")), # max_text_width(genes)),
                      HGVS.c = anno_text(HGVS.c,
                                         location = 0.5, just = "center",
                                         gp = gpar(fontsize = 8, fill = "#D3F7CD", col = "black", border = "white"),
                                         width = unit(5, "cm")), # max_text_width(HGVS.c)),
                      HGVS.p = anno_text(HGVS.p,
                                         location = 0.5, just = "center",
                                         gp = gpar(fontsize = 8, fill = "#F8E5A6", col = "black", border = "white"),
                                         width = unit(5, "cm")), # max_text_width(HGVS.p))
                      pos = anno_text(annot_data$POS,
                                         location = 0.5, just = "center",
                                         gp = gpar(fontsize = 11, fill = "white", col = "black", border = "white"),
                                         width = unit(1.4, "cm"))
) 

data_all <- matrix(0, nrow = nrow(annot_data), ncol = nrow(timepoints))
row.names(data_all) <- annot_data$ID
colnames(data_all) <- timepoints$samples

for (k in timepoints$samples){
  
  temp <- all_variants[which(all_variants$sample == k), c("ID", "AF") ]
  
  data_all[temp$ID, k] <- temp$AF
  
}

colnames(data_all) <- timepoints$heat_labels

if (threshold == 0){
  # save in a csv before 
  final_table <- cbind(data_all, annot_data[, c("POS", "Gene_Name", "HGVS.c", "HGVS.p" ,"REF", "ALT")])
  write.table(final_table, "all variants table.csv", row.names = F, sep = "\t")
}

graph2 <-Heatmap(data_all, name = "frequency", 
                 col = col_fun,
                 heatmap_width = unit(38, "cm"), heatmap_height = unit(nrow(annot_data)/2, "cm"),
                 border = T,
                 #border_gp = gpar(col = "grey"),
                 rect_gp = gpar(col = "white", lwd = 1),
                 cluster_rows = FALSE,cluster_columns  = FALSE,
                 column_names_gp = gpar(fontsize = 13),
                 #row_names_gp = gpar(fontsize = 14),
                 column_names_rot = 50,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", data_all[i,j]), x, y, gp = gpar(fontsize = 12))
                 },
                 heatmap_legend_param = list(
                   labels_gp = gpar(fontsize = 10),
                   title_gp =gpar(fontsize = 13) ,
                   legend_height = unit(2.5, "cm"),
                   legend_width = unit(4, "cm"),
                   at = c(0, threshold,1)
                 ),
                 right_annotation = extra,
                 row_split = factor(genes, levels = unique(genes)),
                 show_row_names = F,
                 #row_order = order(as.numeric(annot_data$POS), decreasing = F),
                 row_title_rot = 0,
                 cluster_row_slices = F,
)


save_image(print(graph2),paste0(plots_folder, "/all variants_", threshold, ".png")
           , width = 20, height = 40, res = 100)



# Specific mutations ----------------------

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

#row.names(heat_matrix) <- new_annotations$mutation
  
extra = rowAnnotation(variants = anno_text(new_annotations$annot, location = 0.5, just = "center",
                                             gp = gpar(fontsize = 11, fill = "#DCE5F3", col = "black", border = "white"),
                                             max_text_width(new_annotations$annot) *1.3) ,
                        Gene_name = anno_text(new_annotations$genes,
                                              location = 0.5, just = "center",
                                              gp = gpar(fontsize = 11, fill = "#F5CDF7", col = "black", border = "white"),
                                              max_text_width(new_annotations$genes) *1.3),
                        HGVS.c = anno_text(new_annotations$HGVS.c,
                                           location = 0.5, just = "center",
                                           gp = gpar(fontsize = 11, fill = "#D3F7CD", col = "black", border = "white"),
                                           max_text_width(new_annotations$HGVS.c) *1.3),
                        HGVS.p = anno_text(new_annotations$HGVS.p,
                                           location = 0.5, just = "center",
                                           gp = gpar(fontsize = 11, fill = "#F8E5A6", col = "black", border = "white"),
                                           max_text_width(new_annotations$HGVS.p) *1.3),
                        pos = anno_text(new_annotations$POS,
                                        location = 0.5, just = "center",
                                        gp = gpar(fontsize = 11, fill = "#ffe0c0", col = "black", border = "white"),
                                        width = max_text_width(new_annotations$POS)),
                        label = anno_text(new_annotations$mutation,
                                        location = 0.5, just = "center",
                                        gp = gpar(fontsize = 11, fill = "white", col = "black", border = "white"),
                                        width = max_text_width(new_annotations$mutation)))


graph3 <-Heatmap(heat_matrix, name = "frequency", 
                   col = col_fun,
                   heatmap_width = unit(31, "cm"), heatmap_height = unit(nrow(heat_matrix) +3, "cm"),
                   border = T,
                   #border_gp = gpar(col = "grey"),
                   rect_gp = gpar(col = "white", lwd = 1),
                   cluster_rows = FALSE,cluster_columns  = FALSE,
                   column_names_gp = gpar(fontsize = 13),
                   #row_names_gp = gpar(fontsize = 14),
                   column_names_rot = 50,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.2f", heat_matrix[i,j]), x, y, gp = gpar(fontsize = 12))
                   },
                   heatmap_legend_param = list(
                     labels_gp = gpar(fontsize = 10),
                     title_gp =gpar(fontsize = 13) ,
                     legend_height = unit(2.5, "cm"),
                     legend_width = unit(4, "cm"),
                     at = c(0, 1)
                   ),
                   right_annotation = extra,
                   #row_split = factor(annotations$genes, levels = unique(annotations$genes)),
                   show_row_names = F,
                   #row_order = order(row.names(heat_matrix), decreasing = T),
                   row_title_rot = 0,
                   cluster_row_slices = F,
  )
  
specific_fp <- str_remove(specific_fp, ".tsv")
specific_fp <- str_remove(specific_fp, "outbreak_info/outbreakinfo_mutation_report_data_")

save_image(print(graph3),paste0(plots_folder, "/",specific_fp,".png"), width = 15, height = 12, res = 100)

# Create dotplot ------------  
if (str_detect(specific_fp, "b.1.1.7")){
  
  create_dotplot(heat_matrix,
                 new_annotations,
                 use_nas = FALSE,
                 plots_folder,
                 specific_fp)
  
}  





