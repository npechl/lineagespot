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

# Read sewage data ---------------------

sewage_file <- "new files/sewage_VCFList.csv"
input_data <- read.csv(sewage_file)
rm(sewage_file)

# Read meta-data file ------------

meta_file <- "new files/meta-data.xlsx"
meta_data <- read.xlsx(meta_file, sheetIndex = 1)

timepoints <- data.table(samples = meta_data$Sample_ID,
                         heat_labels = meta_data$Time.length_.sampling.)

rm(meta_data, meta_file)

# Read VOC/VOI data ---------------

mutations_source <- "outbreak_info" #"pangolin"   #or  "veo"   
specific_fp <- "outbreak_info/"

all_files <- list.files(specific_fp)
all_files <- paste0(specific_fp, all_files)

# Read amino acid abbreviations --------------

amino_acids <- "Amino_acid_Abbreviations.txt"
aa_abbreviations <- read.table(amino_acids, header = T, sep = ",")
rm(amino_acids)

# Create outputs folder -------------

plots_folder <- "new sars/samples_24_5"
dir.create(plots_folder)

# Set threshold for all variants heatmap ----------
threshold <- 0.5

# Create heatmaps ------------------

input_data <- input_data[which(input_data$ID != ""), ]
col_fun = colorRamp2(c( 0,1), c("#FFF6FF", "#BC3C29FF"))

# Common variants --------------------------

#common variants
common_variants <- dplyr::count(input_data, ID)
common_variants <- common_variants[which(common_variants$n == nrow(timepoints)), ]
common_variants <- input_data[which(input_data$ID %in% common_variants$ID ), ]
common_variants$POS <- as.numeric(common_variants$POS)
common_variants <- common_variants[order(common_variants$POS, decreasing = F), ]

if (nrow(common_variants) > 0) {
  
  annot_data <- unique(common_variants[, c("POS" ,"ID", "REF", "ALT", "Gene_Name", "HGVS.c", "HGVS.p")])
  rm(common_variants)
  
  annot_data$annot <- paste0(annot_data$REF, " -> ", annot_data$ALT)
  
  extra = rowAnnotation(variants = anno_text(annot_data$annot, location = 0.5, just = "center",
                                             gp = gpar(fontsize = 11, fill = "#DCE5F3", col = "black", border = "white"),
                                             width = max_text_width(annot_data$annot)*1.1),
                        Gene_name = anno_text( annot_data$Gene_Name,
                                              location = 0.5, just = "center",
                                              gp = gpar(fontsize = 11, fill = "#F5CDF7", col = "black", border = "white"),
                                              width = max_text_width( annot_data$Gene_Name)),
                        HGVS.c = anno_text(annot_data$HGVS.c,
                                           location = 0.5, just = "center",
                                           gp = gpar(fontsize = 11, fill = "#D3F7CD", col = "black", border = "white"),
                                           width = max_text_width(annot_data$HGVS.c)),
                        HGVS.p = anno_text(annot_data$HGVS.p,
                                           location = 0.5, just = "center",
                                           gp = gpar(fontsize = 11, fill = "#F8E5A6", col = "black", border = "white"),
                                           width = max_text_width(annot_data$HGVS.p))
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
                   heatmap_width = unit(18 +nrow(timepoints), "cm"), 
                   heatmap_height = unit(3+nrow(annot_data), "cm"),
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
                   row_split = factor(annot_data$Gene_Name, levels = unique(annot_data$Gene_Name)),
                   row_order = order(as.numeric(row.names(data_common)), decreasing = F),row_title_rot = 0,
                   cluster_row_slices = F,
  )
  
  
  save_image(print(graph1),
             paste0(plots_folder, "/common variants.png"), 
             width = 10, height = 6, res = 100)
  
  output_table <- cbind(data_common, annot_data[, c("annot", "Gene_Name", "HGVS.c", "HGVS.p")])
  
  write.csv(output_table, 
            paste0(plots_folder, "/common_variants.csv"),
            row.names = FALSE)
  
  rm(annot_data, extra, data_common, k, graph1)
}


# All variants --------------------------

#all variants

all_variants <- input_data[which(input_data$AF > threshold), ]
all_variants$POS <- as.numeric(all_variants$POS)
all_variants <-all_variants[order(all_variants$POS, decreasing = F), ]

annot_data <- unique(all_variants[, c("POS" ,"ID", "REF", "ALT", 
                                      "Gene_Name", "HGVS.c", "HGVS.p")])
annot_data$annot <- paste0(annot_data$REF, " -> ", annot_data$ALT)

who <- which(nchar(annot_data$annot) > 20)
annot_data$annot[who] <- paste0(str_sub(annot_data$REF[who],1,10),"...", " -> ",
                                str_sub(annot_data$ALT[who],1,10),"...")

who <- which(nchar(annot_data$HGVS.c) > 25)
annot_data$HGVS.c[who] <- paste0(str_sub(annot_data$HGVS.c[who],1,25),"...")

who <- which(nchar(annot_data$HGVS.p) > 25)
annot_data$HGVS.p[who] <- paste0(str_sub(annot_data$HGVS.p[who],1,25),"...")

rm(who)

extra = rowAnnotation(variants = anno_text(annot_data$annot, location = 0.5, just = "center",
                                           gp = gpar(fontsize = 8, fill = "#DCE5F3", col = "black", border = "white"),
                                           width = unit(5, "cm")), # max_text_width(annot)),
                      Gene_name = anno_text(annot_data$Gene_Name,
                                            location = 0.5, just = "center",
                                            gp = gpar(fontsize = 8, fill = "#F5CDF7", col = "black", border = "white"),
                                            width = unit(2, "cm")), # max_text_width(genes)),
                      HGVS.c = anno_text(annot_data$HGVS.c,
                                         location = 0.5, just = "center",
                                         gp = gpar(fontsize = 8, fill = "#D3F7CD", col = "black", border = "white"),
                                         width = unit(5, "cm")), # max_text_width(HGVS.c)),
                      HGVS.p = anno_text(annot_data$HGVS.p,
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

rm(k, temp)
colnames(data_all) <- timepoints$heat_labels

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
                 row_split = factor(annot_data$Gene_Name, levels = unique(annot_data$Gene_Name)),
                 show_row_names = F,
                 #row_order = order(as.numeric(annot_data$POS), decreasing = F),
                 row_title_rot = 0,
                 cluster_row_slices = F,
)


save_image(print(graph2), 
           paste0(plots_folder, "/all variants_", threshold, ".png"),
           width = 17, height = 30, res = 100)

final_table <- cbind(data_all, annot_data[, c("POS", "annot", "Gene_Name", "HGVS.c", "HGVS.p")])
write.csv(final_table, 
          paste0(plots_folder, "/all variants table", threshold, "csv"), 
          row.names = F)

rm(graph2, data_all, extra, annot_data, all_variants, final_table, threshold)

# Specific mutations ------------------

plots_folder <- paste0(plots_folder, "/", mutations_source)
dir.create(plots_folder)

for (one_file in all_files) {
  
  if (mutations_source == "outbreak_info"){
    
    specific_lineage <- read.table(one_file, header = T, sep = "\t")
    
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
    
    rm(others)
    
    if (nrow(deletions) > 0){
      
      list_del <- check_deletions(input_data,
                                  deletions,
                                  timepoints,
                                  aa_abbreviations)
      
      heat_matrix <- rbind(list_del$data_specific, list_aa$data_specific)
      
      del_annot <- list_del$annotations
      annotations <- rbind(del_annot[, -c("match")], list_aa$annotations)
      row.names(annotations) <- c(del_annot$match, row.names(list_aa$annotations))
      
      rm(deletions, list_del, list_aa, del_annot)
      
    } else {
      
      heat_matrix <- list_aa$data_specific
      annotations <- list_aa$annotations
      
      rm(list_aa)
    }
    
    # row.names(annotations) <- annotations$mutation
    # row.names(heat_matrix) <- annotations$mutation
    
    #heat_matrix <- heat_matrix[specific_lineage$mutation, ]
    #annotations <- annotations[c(row.names(heat_matrix)), ]
    #annotations <- annotations[match(specific_lineage$mutations$V1, row.names(annotations)),]
    annotations$names <- row.names(annotations)
    annotations$names <- reorder.factor(annotations$names, new.order=row.names(heat_matrix))
    new_annotations <- annotations %>%
      arrange(names)
    
    rm(annotations)
    
    new_annotations[is.na(new_annotations)] <- " "
    new_annotations[new_annotations == ""] <- " "
    
    one_file <- str_remove(one_file, ".tsv")
    one_file <- unlist(str_split(one_file, "/"))
    one_file <- one_file[length(one_file)]
    one_file <- str_remove(one_file, "outbreakinfo_mutation_report_data_")  
    
  } else if (mutations_source == "pangolin"){
    
    specific_lineage <- read.table(one_file, header = F, sep = "\t")
    
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
    
    rm(others)
    
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
      
      rm(deletions, list_aa, list_del)
      
    } else {
      
      heat_matrix <- list_aa$data_specific
      annotations <- list_aa$annotations
      rm(list_aa)
    }
    
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
    
    rm(annotations)
    
    new_annotations[is.na(new_annotations)] <- " "
    new_annotations[new_annotations == ""] <- " "
    
    #row.names(heat_matrix) <- new_annotations$mutation
    one_file <- str_remove(one_file, ".txt")
    one_file <- unlist(str_split(one_file, "/"))
    one_file <- one_file[length(one_file)]
    
    new_annotations$mutation <- new_annotations$label
    
  } else if (mutations_source == "veo"){
    
    specific_lineage <- read.table(one_file, header = F, sep = "\t")
    
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
    
    rm(others)
    
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
      
      rm(deletions, list_del, list_aa)
      
    } else {
      
      heat_matrix <- list_aa$data_specific
      annotations <- list_aa$annotations
      rm(list_aa)
    }
    
    row.names(heat_matrix) <- row.names(annotations)
    heat_matrix <- heat_matrix[specific_lineage$label, ]
    
    annotations$label <- reorder.factor(annotations$label, new.order=row.names(heat_matrix))
    new_annotations <- annotations %>%
      arrange(label)
    
    rm(annotations)
    
    new_annotations[is.na(new_annotations)] <- " "
    new_annotations[new_annotations == ""] <- " "
    
    one_file <- str_remove(one_file, ".txt")
    one_file <- unlist(str_split(one_file, "/"))
    one_file <- one_file[length(one_file)]
    
    new_annotations$mutation <- new_annotations$label
    
  }


    
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
  
  heat_width <-  max_text_width(new_annotations$annot) *1.3 +
    max_text_width(new_annotations$genes) *1.3 +
    max_text_width(new_annotations$HGVS.c) *1.3 +
    max_text_width(new_annotations$HGVS.p) *1.3 +
    max_text_width(new_annotations$POS) +
    max_text_width(new_annotations$mutation) 
  #  1.5 * nrow(timepoints) + 6
  
  heat_heigth <- nrow(heat_matrix) + 6
   # heat_heigth <- 27
  graph3 <-Heatmap(heat_matrix, name = "frequency", 
                     col = col_fun,
                     heatmap_width =  heat_width + unit( 1.1 * nrow(timepoints) + 6, "cm"), 
                     heatmap_height = unit(heat_heigth, "cm"),
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
  
  save_as_pdf(print(graph3),
             paste0(plots_folder, "/F9_A.pdf"),
             width = 1.5 * nrow(timepoints) + 10,
             height = heat_heigth/2.5)
  
  save_image(print(graph3),
             paste0(plots_folder, "/",one_file,".png"),
             width = 1.5 * nrow(timepoints) + 10,
             height = heat_heigth/2.5, res = 100)
  
  output_table <- cbind(heat_matrix, 
                        new_annotations[, c("annot", "genes", "HGVS.c", "HGVS.p", "POS", "mutation")])
  
  write.csv(output_table, 
            paste0(plots_folder, "/",one_file,".csv"), 
            row.names = FALSE)
  
  rm(graph3, extra,  
     output_table,  i, specific_lineage)
  
   # Create dotplots ------------  
   
   create_dotplot(heat_matrix,
                  new_annotations,
                  use_nas = TRUE,
                  plots_folder,
                  one_file)
   
   create_dotplot(heat_matrix,
                  new_annotations,
                  use_nas = FALSE,
                  plots_folder,
                  one_file)
     
    
   # Create boxplot ------------
   
   create_boxplot(heat_matrix,
                  new_annotations,
                  plots_folder,
                  one_file)
   
  rm(heat_matrix, new_annotations)

}
