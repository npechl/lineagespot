
rm(list = ls())

# load libraries ----------------------------------

library(data.table)
library(stringr)

# read input data ------------------------

meta = readxl::read_excel("meta-data.xlsx", sheet = 1)

sewage_data = fread("sewage_VCFList.csv")

# remove "p." character and correct ORF ----------------

sewage_data$HGVS.p = str_remove(sewage_data$HGVS.p, "p\\.") 


# clean AA variants ---------------------------------------

aa_variants = str_split(sewage_data$HGVS.p, "[0-9]|_")
codon_num = str_split(sewage_data$HGVS.p, "[A-Z]|[a-z]|\\_|\\*|\\?")

aa_variants = lapply(aa_variants, function(x) {
  
  return( x[x != ""] )
  
})

codon_num = lapply(codon_num, function(x) {
  
  return( x[x != ""] )
  
})


# merge overall tables ---------------------------------------

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


# add AA variants with Nt when there is not any -------------------------

sewage_data[which(sewage_data$HGVS.p == ""), ]$HGVS.p = sewage_data[which(sewage_data$HGVS.p == ""), ]$ID


# clean AA abbreviation --------------------------------

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


voc_data = sewage_data

voc_data = unique(voc_data)



# collapse table ------------------------

voc_data = voc_data[, c("POS", "gt_DP", "gt_AD_alt", 
                        "Gene_Name", "HGVS.p", "sample"), 
                    with = FALSE]



voc_data = voc_data[, .(DP = unique(gt_DP),
                        AD = sum(gt_AD_alt)), 
                    by = .(POS, Gene_Name, HGVS.p, sample)]

voc_positions = unique(voc_data[, c("POS", "HGVS.p", "Gene_Name"), 
                                with = FALSE])

voc_data = voc_data[, .(DP = sum(DP),
                        AD = sum(AD)), 
                    by = .(Gene_Name, HGVS.p, sample)]


voc_data$AF = voc_data$AD / voc_data$DP

voc_data$barcode = paste0(voc_data$Gene_Name, ":", voc_data$HGVS.p)


who = match(voc_data$sample, meta$Sample_ID)

voc_data$time_point = meta[who, ]$`Time length_(sampling)`

voc_data$time_point = factor(voc_data$time_point, 
                             levels = meta$`Time length_(sampling)`)

rm(list = setdiff(ls(), c("voc_data", "meta", "voc_positions")))





  
# Figure 2.A -----------------------------------------------

plot_data = voc_data[, c("sample", "AF", "barcode"), with = FALSE]

plot_data = plot_data %>% tidyr::spread(key = "sample", value = "AF", fill = 0)

barcodes = plot_data$barcode

plot_data = as.matrix(plot_data[, meta$Sample_ID, with = FALSE])

colnames(plot_data) = meta$`Time length_(sampling)`
row.names(plot_data) = barcodes


dend = dist(plot_data, method = "euclidean")

dend = hclust(dend, method = "ward.D")

den = as.dendrogram(dend)


library(ggtree)

clus = cutree(dend, 2)

g = split(names(clus), clus)

p = ggtree(dend)

clades = sapply(g, function(n) MRCA(p, n))

p = groupClade(p, clades, group_name = 'subtree') + aes(color = subtree)

p + layout_dendrogram() +
    
  scale_color_manual(values = c("blue", "gray50", "red")) + 
    
  theme_dendrogram(plot.margin = margin(6,6,80,6)) +
    
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())





# Figure 2.B ----------------------------------------------

library(ComplexHeatmap)

split_tree = cutree(dend, k = 5)

table(split_tree)

who = names(split_tree[which(split_tree %in% c(5))])

temp = plot_data[who, ]

ht1 = Heatmap(matrix = temp,
              name = "Variant evolution",
              
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D2",
              
              col = c("#91D1C299", "yellow2"),
              
              cluster_columns = FALSE,
              
              show_row_dend = TRUE,
              
              show_row_names = TRUE,
              
              row_title = NULL,
              
              border = FALSE,
              
              row_dend_width = unit(20, "mm"),
              
              show_heatmap_legend = TRUE,
              
              heatmap_legend_param = list(legend_height = unit(3, "cm")),
              
              rect_gp = gpar(col = "white", lwd = 1.2))


draw(ht1)

rm(list = setdiff(ls(), c("voc_data", "meta", "voc_positions")))

# Figure 2.C --------------------------------

plot_data = voc_data[, .N, by = .(Gene_Name, sample)]

plot_data = plot_data[which(plot_data$Gene_Name != "ORF7b&ORF8"), ]

gene_meta = readxl::read_excel("gene-metadata.xlsx", sheet = 1)


who = match(plot_data$Gene_Name, gene_meta$Gene_name)

plot_data$N = round(plot_data$N / gene_meta[who, ]$`Gene_length (nt)`, 
                    digits = 4)


who = match(plot_data$sample, meta$Sample_ID)
plot_data$time_point = meta[who, ]$`Time length_(sampling)`

plot_data$time_point = factor(plot_data$time_point, 
                              levels = meta$`Time length_(sampling)`)

who = as.numeric(plot_data$time_point)

plot_data$highlight_point = "0.8"

plot_data[which(who == 9), ]$highlight_point = "1"

plot_data$highlight_line = "0.8"

plot_data[which(who %in% c(8, 9)), ]$highlight_line = "1"

gc()

library(ggplot2)

ggplot(plot_data, 
       aes(x = time_point, y = N, group = Gene_Name)) + 
  
  geom_vline(xintercept = "19-25 March 2021", linetype="dotted",
             color = "gray90", size=0.9) +

  geom_vline(xintercept = "2-8 April 2021", linetype="dotted",
             color = "gray90", size=0.9) +
  
  geom_line(aes(alpha = highlight_line),
            size = 1) +
  
  
  
  # scale_color_manual(values = c("gray70", "gray10")) +
  
  geom_point(size = 2, aes(alpha = highlight_point)) +
  
  scale_alpha_manual(values = c(0.3, 1)) +
    
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, 0.15)) +
  
  theme_minimal() +
  
  theme(# panel.grid = element_blank(),
        # panel.border = element_rect(fill = NA),
    
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
    
        legend.position = "none",
        
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                                   size = 10, face = "bold"),
        
        axis.title.x = element_blank(),
        
        axis.line.x = element_line(),
        
        axis.ticks = element_line(),
        # axis.ticks.y = element_line(),
        
        axis.text.y = element_text(size = 8, face = "bold"),
        
        axis.title.y = element_blank(),
        
        axis.line.y = element_blank(),
        
        strip.placement = "outside",
        
        strip.text = element_text(face = "bold")
        ) + 
  
  facet_wrap(~Gene_Name, ncol = 1, strip.position = "left", scales = "free_y") +
  
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  
  labs(title = "Normalized number of mutations per gene", linetype = "Gene")



citation("data.table")
citation("stringr")
citation("ggtree")
citation("ComplexHeatmap")
citation("ggplot2")



