
rm(list = ls())

gc()

citation("data.table")
citation("stringr")
citation("ComplexHeatmap")

# load libraries ------------------------------

library(data.table)
library(stringr)

# read input data -----------------------------------

meta = readxl::read_excel("meta-data.xlsx", sheet = 1)

sewage_data = fread("sewage_VCFList.csv")

# remove "p." character and correct ORF --------------------

sewage_data$HGVS.p = str_remove(sewage_data$HGVS.p, "p\\.") 


# clean AA variants -------------------------

aa_variants = str_split(sewage_data$HGVS.p, "[0-9]|_")
codon_num = str_split(sewage_data$HGVS.p, "[A-Z]|[a-z]|\\_|\\*|\\?")

aa_variants = lapply(aa_variants, function(x) {
  
  return( x[x != ""] )
  
})

codon_num = lapply(codon_num, function(x) {
  
  return( x[x != ""] )
  
})


# merge overall tables --------------------------

aa_split_list = lapply(seq(1:nrow(sewage_data)), function(i) {
  
  out = data.table(ref_aa = NA,
                   alt_aa = NA,
                   info = NA,
                   codon_num = NA,
                   codon_end = NA)
  
  
  if(length(aa_variants[[i]]) == 2) {
    
    out$ref_aa = aa_variants[[i]][1]
    out$alt_aa = aa_variants[[i]][2]
    
    out$codon_num = codon_num[[i]][1]
    
  } else if(length(aa_variants[[i]] == 3)) {

    out$ref_aa = aa_variants[[i]][1]
    out$alt_aa = aa_variants[[i]][2]
    
    out$type == aa_variants[[i]][3]
    
    out$codon_num = codon_num[[i]][1]
    out$codon_end = codon_num[[i]][2]
    
  } else {
    
    
  }
  
  return(out)
  
})

aa_split_list = rbindlist(aa_split_list)

sewage_data = cbind(sewage_data, aa_split_list)

sewage_data$change_length_nt = abs(str_length(sewage_data$ALT) - 
                                       str_length(sewage_data$REF))

sewage_data$codon_num = as.numeric(sewage_data$codon_num)
sewage_data$codon_end = as.numeric(sewage_data$codon_end)


# add AA variants with Nt when there is not any ---------------

sewage_data[which(sewage_data$HGVS.p == ""), ]$HGVS.p = sewage_data[which(sewage_data$HGVS.p == ""), ]$ID


# clean AA abbreviation --------------------------

AA_abreviations = fread("Amino_acid_Abbreviations.txt")

for(i in 1:nrow(AA_abreviations)) {
  
  sewage_data$ref_aa = str_replace_all(sewage_data$ref_aa, 
                                       AA_abreviations[i,]$Three_Letter,
                                       AA_abreviations[i,]$One_Letter)
  
  sewage_data$alt_aa = str_replace_all(sewage_data$alt_aa, 
                                       AA_abreviations[i,]$Three_Letter,
                                       AA_abreviations[i,]$One_Letter)
  
  sewage_data$HGVS.p = str_replace_all(sewage_data$HGVS.p,
                                       AA_abreviations[i,]$Three_Letter,
                                       AA_abreviations[i,]$One_Letter)
  
}


# read reference -----------------------------------

ref = fread("outbreakinfo_mutation_report_data_b.1.1.7.tsv")

ref[which(ref$gene == "ORF1a"), ] $gene = "ORF1ab"
ref[which(ref$gene == "ORF1b"), ] $gene = "ORF1ab"

ref = ref[order(ref$gene), ]


ref$barcode = paste0(ref$gene, ":", ref$ref_aa, ref$codon_num, ref$alt_aa)

ref[which(ref$type == "deletion"), ]$barcode = ref[which(ref$type == "deletion"), ]$ref_aa

ref[which(ref$change_length_nt == "None"), ]$change_length_nt = "0"

ref$change_length_nt = as.numeric(ref$change_length_nt)


# B.1.1.7 hits -----------------------------

ref$alt_aa = str_replace(ref$alt_aa, "\\*", "\\\\*")
sewage_data$alt_aa = str_replace(sewage_data$alt_aa, "\\*", "\\\\*")

voc_data = list()

for(i in 1:nrow(ref)) {
  
  
  if(ref[i,]$type == "substitution") {
    
    temp = sewage_data[which(sewage_data$Gene_Name == ref[i,]$gene), ]
    
    temp = temp[which(str_detect(temp$ref_aa, 
                                 ref[i,]$ref_aa)), ]
    
    temp = temp[which(str_detect(temp$alt_aa, 
                                 ref[i,]$alt_aa)), ]
    
    temp = temp[which(temp$codon_num == ref[i,]$codon_num | 
                          temp$codon_num == ref[i,]$codon_num + 1 | 
                          temp$codon_num == ref[i,]$codon_num - 1), ]
    
    temp = temp[which(temp$change_length_nt == ref[i,]$change_length_nt), ]
    
  } else {
    
    temp = sewage_data[which(sewage_data$Gene_Name == ref[i,]$gene), ]
    temp = temp[which(str_detect(temp$HGVS.p, "del")), ]
    
    temp = temp[which(temp$codon_num == ref[i,]$codon_num), ]
    temp = temp[which(temp$codon_end == ref[i,]$codon_end), ]
    
    temp = temp[which(temp$change_length_nt == ref[i,]$change_length_nt), ]
    
    
  }
  
  
  voc_data[[ref[i,]$barcode]] = temp
  

  
}

not_overlapping_variants = names(which(lapply(voc_data, nrow) == 0))


voc_data = rbindlist(voc_data)

voc_data = unique(voc_data)



# collapse table -------------------------------

voc_data = voc_data[, c("POS", "gt_DP", "gt_AD_alt",
                        "Gene_Name", "HGVS.p", "sample"), 
                    with = FALSE]



voc_data = voc_data[, .(DP = unique(gt_DP),
                        AD = sum(gt_AD_alt)), 
                    by = .(POS, Gene_Name, HGVS.p, sample)]

voc_positions = unique(voc_data[, c("POS", "HGVS.p", "Gene_Name"), with = FALSE])

voc_data = voc_data[, .(DP = sum(DP),
                        AD = sum(AD)), 
                    by = .(Gene_Name, HGVS.p, sample)]


voc_data$AF = voc_data$AD / voc_data$DP


# add non overlapping rules ---------------------


not_overlapping_variants = str_split(not_overlapping_variants, "\\:", simplify = TRUE)

voc_data = rbind(voc_data, 
                 data.table(Gene_Name = not_overlapping_variants[,1],
                            HGVS.p = not_overlapping_variants[,2],
                            sample = unique(voc_data$sample),
                            DP = 0,
                            AD = 0,
                            AF = 0))

voc_positions = rbind(voc_positions,
                      data.table(POS = c(14368, 21991),
                                 HGVS.p = not_overlapping_variants[,2],
                                 Gene_Name = not_overlapping_variants[,1]))

voc_positions = unique(voc_positions)


rm(list = setdiff(ls(), c("voc_data", "meta", "voc_positions")))


# Plotting part -------------------------------


voc_positions$barcode = paste0(voc_positions$Gene_Name, ":", voc_positions$HGVS.p)

voc_positions = voc_positions[, c("POS", "barcode"), with = FALSE]


coverage_data = fread("coverage-all-sewage-samples.csv")

voc_coverage = coverage_data[voc_positions$POS, ]

voc_coverage$barcode = voc_positions$barcode

voc_coverage = voc_coverage[, lapply(.SD, sum), by = barcode]

rm(coverage_data)


voc_data$barcode = paste0(voc_data$Gene_Name, ":", voc_data$HGVS.p)

plot_matrix = voc_data[, c("sample", "barcode", "AF"), with = FALSE]

plot_matrix = plot_matrix %>% tidyr::spread(key = "sample", value = "AF", fill = 0)


plot_matrix = plot_matrix[, c("barcode", meta$Sample_ID), with = FALSE]
voc_coverage = voc_coverage[, c("barcode", meta$Sample_ID), with = FALSE]

plot_matrix = plot_matrix[order(plot_matrix$barcode), ]
voc_coverage = voc_coverage[order(voc_coverage$barcode), ]


temp = plot_matrix$barcode

plot_matrix = as.matrix(plot_matrix[,2:ncol(plot_matrix)])

row.names(plot_matrix) = temp
colnames(plot_matrix) = meta$`Time length_(sampling)`

temp = voc_coverage$barcode

voc_coverage = as.matrix(voc_coverage[,2:ncol(voc_coverage)])

row.names(voc_coverage) = temp
colnames(voc_coverage) = meta$`Time length_(sampling)`

# voc_coverage = log(voc_coverage)

genes = str_split(temp, "\\:", simplify = TRUE)[,1]

# voc_coverage = log(voc_coverage + 1)


library(ComplexHeatmap)

ha = HeatmapAnnotation(`Variants' evolution` = anno_boxplot(plot_matrix),
                       annotation_name_side = "right",
                       annotation_name_rot = 0,
                       annotation_name_offset = unit(2, "mm"))

ha2 = HeatmapAnnotation(`Coverage summary` = anno_density(log10(voc_coverage + 1),
                                                 type = "violin",
                                                 
                                                 axis_param = list(at = log10(c(10, 1000, 100000)),
                                                                   labels = c(10, 1000, 100000))),
                        
                        
                        annotation_name_side = "right",
                        annotation_name_rot = 0,
                        annotation_name_offset = unit(2, "mm"))

ha3 = rowAnnotation(Gene = anno_block(gp = gpar(fill = 2:5),
                                      labels = unique(genes), 
                                      labels_gp = gpar(col = "white", fontsize = 10)))

lgd_list = list( Legend(labels = c("less than 20 reads coverage"),
                        title = "Coverage",
                        type = "grid",
                        # pch = 16,
                        legend_gp = gpar(fill = c("gray20"))))
  


ht2 = Heatmap(matrix = plot_matrix,
              name = "Allele frequency",
              
              cell_fun = function(j, i, x, y, width, height, fill) {
                

                
                # if(voc_coverage[i, j] > 20) {
                #   
                #   grid.polygon( unit.c(x - 0.5 * width, x - 0.5 * width, x + 0.5 * width), 
                #                 unit.c(y - 0.5 * height, y + 0.5 * height, y - 0.5 * height),
                #                 
                #                 gp = gpar(col = "white", fill = "gray90"))
                #   
                # } else {
                #   
                #   grid.polygon( unit.c(x - 0.5 * width, x - 0.5 * width, x + 0.5 * width), 
                #                 unit.c(y - 0.5 * height, y + 0.5 * height, y - 0.5 * height),
                #                 
                #                 gp = gpar(col = "white", fill = "gray20"))
                #   
                # }
                
                grid.text(sprintf("%.1f", plot_matrix[i, j]),
                          x, y,
                          gp = gpar(fontsize = 10))
                
                if(voc_coverage[i, j] < 20) {

                  grid.rect(x, y, width, height, gp = gpar(col = "gray20", # "#E64B35FF",
                                                           fill = "transparent"))

                }
                
              },
              
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D2",
              
              column_title = "Alpha variant (B.1.1.7) evolution",
              
              col = c("#91D1C299", "yellow2"),
        
              cluster_columns = FALSE,
              
              show_row_dend = FALSE,
              
              # row_title = NULL,
              
              row_split = genes,
              
              border = FALSE,
              
              top_annotation = ha,
              
              bottom_annotation = ha2,
              
              # left_annotation = ha3,
              
              # column_names_rot = 45,
              
              show_heatmap_legend = TRUE,
              
              rect_gp = gpar(col = "white", lwd = 1.2),
              
              heatmap_legend_param = list(legend_width = unit(3, "cm"),
                                          direction = "horizontal")
              )

draw(ht2, 
     annotation_legend_list = lgd_list, 
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legend = TRUE)


