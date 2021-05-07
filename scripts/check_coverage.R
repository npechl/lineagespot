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

# Read coverage data ----------------------

coverage_file <- "coverage-all-sewage-samples.csv"
coverage_data <- read.csv(coverage_file)
rm(coverage_file)

# Read sewage data ---------------

sewage_file <- "new files/sars-sewage.csv"
input_data <- read.csv(sewage_file)
rm(sewage_file)

# Read meta-data file ------------

meta_file <- "new files/meta-data.xlsx"
meta_data <- read.xlsx(meta_file, sheetIndex = 1)

timepoints <- data.table(samples = meta_data$Sample_ID,
                         heat_labels = meta_data$Time.length_.sampling.)

rm(meta_data, meta_file)

# Read VOC and VOI files ------------

mutations_source <-  "pangolin"   # "veo" #or  "outbreak_info" 
specific_fp <- "pangolin data/b.1.351_pangolin.txt"

# Read amino acid abbreviations file -------------

amino_acids <- "Amino_acid_Abbreviations.txt"
aa_abbreviations <- read.table(amino_acids, header = T, sep = ",")
rm(amino_acids)

# Create output folder ---------------

plots_folder <- "new sars/new samples/coverage"
dir.create(plots_folder)

# Read gene coordinates -------------

gene_coordinates <- "new files/NC_045512.2_annot.gff3"
file_genes <- fread(gene_coordinates, header = FALSE, sep = "\t")
rm(gene_coordinates)

# Create heatmaps ------------------------

input_data <- input_data[which(input_data$ID != ""), ]
col_fun = colorRamp2(c( 0,1), c("#FFF6FF", "red"))

if (mutations_source == "outbreak_info"){
  
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
  
  annotations$names <- row.names(annotations)
  annotations$names <- reorder.factor(annotations$names, new.order=row.names(heat_matrix))
  new_annotations <- annotations %>%
    arrange(names)
  
  rm(annotations)
  
  new_annotations[is.na(new_annotations)] <- " "
  new_annotations[new_annotations == ""] <- " "
  
  specific_fp <- str_remove(specific_fp, ".tsv")
  specific_fp <- str_split(specific_fp, "/") %>%
    unlist()
  specific_fp <- specific_fp[length(specific_fp)]
  specific_fp <- str_remove(specific_fp, "outbreakinfo_mutation_report_data_")  
  
} else if (mutations_source == "pangolin"){
  
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
  
  specific_fp <- str_remove(specific_fp, ".txt")
  specific_fp <- unlist(str_split(specific_fp, "/"))
  specific_fp <- specific_fp[length(specific_fp)]
  
  new_annotations$mutation <- new_annotations$label
  
} else if (mutations_source == "veo"){
  
  specific_lineage <- read.table(specific_fp, header = T, sep = "\t")
  
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
  
  specific_fp <- str_remove(specific_fp, ".txt")
  specific_fp <- unlist(str_split(specific_fp, "/"))
  specific_fp <- specific_fp[length(specific_fp)]
  
  new_annotations$mutation <- new_annotations$label
  
}

rm(specific_lineage, input_data)

positions <- new_annotations$HGVS.p
positions <- str_remove_all(positions, paste(aa_abbreviations$Three_Letter, collapse ="|"))
positions <- str_remove_all(positions, "p.|del|\\*| ")

genes <- unlist(str_split(new_annotations$mutation, ":"))

if (mutations_source == "outbreak_info"){
  genes <- genes[c(T, F)]
} else if (mutations_source == "pangolin") {
  genes <- genes[c(F,T,F)]
  genes <- toupper(genes)
  
}

pos_table <- data.table(gene = genes,
                        positions = positions,
                        start = positions,
                        end = positions)

rm(positions, genes)

who <- str_detect(pos_table$positions, "_")
pos_table$start[who] <- unlist(str_split(pos_table$start[who], "_"))[c(T, F)] 
pos_table$end[who] <- unlist(str_split(pos_table$end[who], "_"))[c(F, T)]
rm(who)

pos_table$start <- as.numeric(pos_table$start)
pos_table$end <- as.numeric(pos_table$end)

pos_table$gene <- toupper(pos_table$gene)

gene_annot <- data.table(start = file_genes$V4,
                         end = file_genes$V5,
                         gene = file_genes$V9)

one_row <- data.table(start = gene_annot$start[1],
                      end = gene_annot$end[16],
                      gene = "ORF1AB")

gene_annot <- rbind(one_row, gene_annot[17:nrow(gene_annot), ])
rm(one_row)

gene_annot$gene <- toupper(gene_annot$gene)
gene_annot$start <- as.numeric(gene_annot$start)
gene_annot$end <- as.numeric(gene_annot$end)

pos_table$new_start <- NA
pos_table$new_end <- NA


if (mutations_source == "outbreak_info") {
  pos_table$gene <- str_replace_all(pos_table$gene, "ORF1A|ORF1B", "ORF1AB")
}


for (i in c(1:nrow(pos_table))) {
  
  gene_row <- gene_annot[which(gene_annot$gene == pos_table$gene[i]), ]
  pos_table$new_start[i] <- gene_row$start[1] + ((pos_table$start[i] - 1) * 3) 
  pos_table$new_end[i] <- gene_row$start[1] + ((pos_table$end[i] - 1) * 3) + 2
}

rm(gene_row)

if (mutations_source == "pangolin") {
  
  who_del <- which(str_detect(new_annotations$mutation, "del"))
  pos_table$new_start[who_del] <- as.numeric(pos_table$gene[who_del])
  temp <- unlist(str_split(new_annotations$mutation[who_del], ":"))[c(F,F,T)]
  pos_table$new_end[who_del] <- as.numeric(pos_table$gene[who_del]) + as.numeric(temp)
  new_annotations$names <- new_annotations$mutation
  
}


pos_table$names <- new_annotations$names
who_del <- which(str_detect(new_annotations$names, "del"))
pos_table$new_start[who_del] <- pos_table$new_start[who_del] - 1
pos_table$new_end[who_del] <- pos_table$new_start[who_del]
rm(who_del)

heat_coverage <- matrix(0, nrow = nrow(heat_matrix), ncol = ncol(heat_matrix))

for (i in 1:nrow(pos_table)){
  
  seq_row <- seq(pos_table$new_start[i], pos_table$new_end[i], 1)
  temp <- coverage_data[which(row.names(coverage_data) %in% seq_row), ]
  temp <- colSums(temp)
  heat_coverage[i, ] <- temp
}

rm(seq_row, temp, aa_abbreviations, file_genes, gene_annot, pos_table)

who_aa <- which(str_detect(new_annotations$names, "del", negate = TRUE))
heat_coverage[who_aa, ] <- round(heat_coverage[who_aa, ] / 3)
rm(who_aa)

extra = rowAnnotation( #variants = anno_text(new_annotations$annot, location = 0.5, just = "center",
#                                            gp = gpar(fontsize = 11, fill = "#DCE5F3", col = "black", border = "white"),
#                                            max_text_width(new_annotations$annot) *1.3) ,
#                       Gene_name = anno_text(new_annotations$genes,
#                                             location = 0.5, just = "center",
#                                             gp = gpar(fontsize = 11, fill = "#F5CDF7", col = "black", border = "white"),
#                                             max_text_width(new_annotations$genes) *1.3),
#                       HGVS.c = anno_text(new_annotations$HGVS.c,
#                                          location = 0.5, just = "center",
#                                          gp = gpar(fontsize = 11, fill = "#D3F7CD", col = "black", border = "white"),
#                                          max_text_width(new_annotations$HGVS.c) *1.3),
#                       HGVS.p = anno_text(new_annotations$HGVS.p,
#                                          location = 0.5, just = "center",
#                                          gp = gpar(fontsize = 11, fill = "#F8E5A6", col = "black", border = "white"),
#                                          max_text_width(new_annotations$HGVS.p) *1.3),
#                       pos = anno_text(new_annotations$POS,
#                                       location = 0.5, just = "center",
#                                       gp = gpar(fontsize = 11, fill = "#ffe0c0", col = "black", border = "white"),
#                                       width = max_text_width(new_annotations$POS)),
                      label = anno_text(new_annotations$mutation,
                                        location = 0.5, just = "center",
                                        gp = gpar(fontsize = 11, fill = "white", col = "black", border = "white"),
                                        width = max_text_width(new_annotations$mutation)))


graph3 <-Heatmap(heat_matrix, name = "frequency", 
                 col = col_fun,
                 heatmap_width = unit(23, "cm"), heatmap_height = unit(nrow(heat_matrix) + 11, "cm"),
                 border = T,
                 #border_gp = gpar(col = "grey"),
                 rect_gp = gpar(col = "white", lwd = 1),
                 cluster_rows = FALSE,cluster_columns  = FALSE,
                 column_names_gp = gpar(fontsize = 13),
                 #row_names_gp = gpar(fontsize = 14),
                 column_names_rot = 50,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%d", heat_coverage[i,j]), x, y, gp = gpar(fontsize = 12))
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
 
 
save_image(print(graph3),
           paste0(plots_folder, "/coverage_1",specific_fp,".png"), 
           width = 12, height = 10, res = 300)

colnames(heat_coverage) <- paste0("Coverage_", timepoints$samples)
output_table <- cbind(heat_coverage, heat_matrix, new_annotations[, 1:6])

write.csv(output_table, 
          paste0(plots_folder, "/coverage_",specific_fp,".csv"), 
          row.names = FALSE)

rm(graph3, extra, heat_coverage, heat_matrix, new_annotations, 
   i, output_table, mutations_source)

# Create violin + box plot --------------

# create a dataset
 data <- data.frame(
   name=c( rep("0",nrow(coverage_data)),
           rep("1",nrow(coverage_data)),
           rep("2",nrow(coverage_data)),
           rep("3",nrow(coverage_data)),
           rep("4",nrow(coverage_data)),
           rep("5",nrow(coverage_data)),
           rep("6",nrow(coverage_data)),
           rep("7",nrow(coverage_data)),
           rep("8",nrow(coverage_data)),
           rep("9",nrow(coverage_data)),
           rep("91",nrow(coverage_data))),
   hits=as.vector(unlist(coverage_data))
 )
 
# Violin Plot
 plot <- data %>%
   #left_join(sample_size) %>%
   #mutate(myaxis = paste0(name, "\n", "n=", num)) %>%
   ggplot( aes(x=name, y=hits, fill=name)) +
   geom_violin(width=1.4) +
   geom_boxplot(width=0.1, color="black", alpha=0) +
   scale_fill_viridis(discrete = TRUE) +
   theme_ipsum() +
   theme(
     legend.position="none",
     plot.title = element_text(size=18),
     axis.text.x = element_text(colour = "black", size = 15, angle = 90),
     axis.text.y = element_text(colour = "black", size = 16)
   ) +
   ggtitle("Coverage") +
   xlab("") +
   scale_x_discrete(labels=timepoints$samples, 
                    breaks = c("0","1","2","3","4","5","6", "7", "8", "9", "91")) 
 
 save_image(print(plot),
            paste0(plots_folder, "/coverage_violin.png"), 
            width = 20, height = 10, res = 200)

rm(timepoints, coverage_data, data, plot)
