

rm(list = ls())

# load libraries ---------------------------------------------------------------

library(data.table)
library(stringr)
library(readxl)

library(tidyr)

library(ComplexHeatmap)

# read input data --------------------------------------------------------------

sewage_data = fread("../Thessaloniki/sewage_VCFList.csv")

meta_data = readxl::read_excel("../Thessaloniki/meta-data.xlsx")


# collapse on Amino Acid level -------------------------------------------------

sewage_data[which(sewage_data$HGVS.p == ""), ]$HGVS.p = paste(sewage_data[which(sewage_data$HGVS.p == ""), ]$Gene_Name,
                                                              sewage_data[which(sewage_data$HGVS.p == ""), ]$HGVS.c,
                                                              sep = ";")

sewage_data[which(sewage_data$HGVS.p != ""), ]$HGVS.p = paste(sewage_data[which(sewage_data$HGVS.p != ""), ]$Gene_Name,
                                                              sewage_data[which(sewage_data$HGVS.p != ""), ]$HGVS.p,
                                                              sep = "_")

amino_variants = sewage_data[, .(DP_sum = sum(gt_DP), 
                                 AD_sum = sum(gt_AD_alt)), by = .(sample, HGVS.p)]


amino_variants$AF = amino_variants$AD_sum / amino_variants$DP_sum


# complex heatmap of all AA variants -------------------------------------------

sewage_matrix = amino_variants[,c("HGVS.p", "AF", "sample"), with = FALSE]

sewage_matrix = sewage_matrix %>% spread(key = "sample", value = "AF", fill = 0)

plot_matrix = as.matrix(sewage_matrix[, meta_data$Sample_ID, with = FALSE])

colnames(plot_matrix) = meta_data$`Time length_(sampling)`
row.names(plot_matrix) = sewage_matrix$HGVS.p


dend = dist(plot_matrix, method = "euclidean")

dend = hclust(dend, method = "ward.D")

split_tree = cutree(dend, k = 5)

gene_annotation = rowAnnotation(`Gene name` = str_split(row.names(plot_matrix), "_", simplify = TRUE)[,1],
                                col = list(`Gene name` = c("E" = "#004c4c",
                                                           "M" = "#006666",
                                                           "N" = "#008080",
                                                           "ORF10" = "#66b2b2",
                                                           "ORF1ab" = "blue4",
                                                           "ORF3a" = "#FFFFBF",
                                                           "ORF6" = "#E6F598",
                                                           "ORF7a" = "#ABDDA4",
                                                           "ORF7b" = "#66C2A5",
                                                           "ORF8" = "#3288BD",
                                                           "ORF7b&ORF8" = "black",
                                                           "S" = "gray80")))

ComplexHeatmap::Heatmap(plot_matrix,
                        name = "AA variant frequency",
                        col = c("#0072B5FF", "white", "#BC3C29FF"),
                        rect_gp = gpar(col = "gray12", lty = 0.2, lwd = 0.02),
                        
                        cluster_rows = dend,
                        
                        row_dend_width = unit(3.6, "cm"),
                        row_dend_gp = gpar (15),
                        
                        right_annotation = gene_annotation,
                        
                        row_split = 5,
                        
                        cluster_row_slices = FALSE,
                        row_title_rot = 0,
                        
                        cluster_columns = FALSE,
                        show_row_names = FALSE,
                        column_names_rot = 45)


# complex heatmap of specific clusters variants --------------------------------

dend = dist(plot_matrix, method = "euclidean")

dend = hclust(dend, method = "ward.D")

split_tree = cutree(dend, k = 5)

table(split_tree)

# who = which(as.vector(apply(plot_matrix, 1, var)) >= 0.02)
who = names(split_tree[which(split_tree %in% c(3, 5))])

temp = plot_matrix[who, ]

dend = dist(temp, method = "euclidean")

dend = hclust(dend, method = "ward.D")

gene_annotation = rowAnnotation(`Gene name` = str_split(row.names(temp), "_", simplify = TRUE)[,1],
                                col = list(`Gene name` = c("E" = "#004c4c",
                                                           "M" = "#006666",
                                                           "N" = "#008080",
                                                           "ORF10" = "#66b2b2",
                                                           "ORF1ab" = "blue4",
                                                           "ORF3a" = "#FFFFBF",
                                                           "ORF6" = "#E6F598",
                                                           "ORF7a" = "#ABDDA4",
                                                           "ORF7b" = "#66C2A5",
                                                           "ORF8" = "#3288BD",
                                                           "ORF7b&ORF8" = "black",
                                                           "S" = "gray80")))

ComplexHeatmap::Heatmap(temp,
                        name = "AA variant frequency",
                        col = c("#0072B5FF", "white", "#BC3C29FF"),
                        rect_gp = gpar(col = "gray12", lty = 0.2, lwd = 0.02),
                        
                        cluster_rows = dend,
                        
                        row_dend_width = unit(3.6, "cm"),
                        row_dend_gp = gpar (15),
                        
                        right_annotation = gene_annotation,
                        
                        row_split = 2,
                        
                        cluster_row_slices = FALSE,
                        row_title_rot = 0,
                        
                        cluster_columns = FALSE,
                        show_row_names = FALSE,
                        row_names_side = "right",
                        row_names_gp = gpar(fontsize = 3),
                        column_names_rot = 45)


# gene counts ------------------------------------------------------------------

amino_variants$gene_name = str_split(amino_variants$HGVS.p, "_", simplify = TRUE)[,1]


gene_counts = amino_variants[, .N, by = .(sample, gene_name)]

gene_counts = gene_counts[which(gene_counts$gene_name != "ORF7b&ORF8"), ]

gene_meta = readxl::read_excel("gene-metadata.xlsx", sheet = 1)

who = match(gene_counts$gene_name, gene_meta$Gene_name)

gene_counts$N = 100 * gene_counts$N / gene_meta[who, ]$`Gene_length (nt)`

gene_counts = gene_counts %>% spread(key = "sample", value = "N", fill = 0)


plot_matrix = as.matrix(gene_counts[,meta_data$Sample_ID, with = FALSE])

colnames(plot_matrix) = meta_data$`Time length_(sampling)`
row.names(plot_matrix) = gene_counts$gene_name

ComplexHeatmap::Heatmap(plot_matrix,
                        name = "Gene variants",
                        col = c("#0072B5FF", "white", "#BC3C29FF"),
                        rect_gp = gpar(col = "gray40", lty = 1, lwd = 1),
                        
                        clustering_distance_rows = "euclidean",
                        clustering_method_rows = "ward.D",
                        
                        row_dend_width = unit(3.6, "cm"),
                        row_dend_gp = gpar (15),
                        
                        row_title_rot = 0,
                        
                        cluster_columns = FALSE,
                        show_row_names = TRUE,
                        column_names_rot = 45)









