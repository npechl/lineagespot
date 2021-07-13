

rm(list = ls())

# load libraries ---------------------------------------------------------------

source("libraries.R")
source("helpers.R")

library(ggplot2)
library(ggsci)

# Read sample list -----------------

samplelist = readxl::read_excel("../Replicate-samples/meta-data.xlsx")

freq_values = c("01", "10", "20", "30", "40", "50")

total_vcf_list = list()

for(f in freq_values) {
  
  vcf_list = list()
  
  for(i in samplelist[[1]]) {
    
    vcf_list[[i]] = read.vcfR(paste0("../Replicate-samples/", i, "_freebayes_", f, "_ann.vcf"), 
                              verbose = FALSE)
    
  }
  
  rm(i)
  
  
  # Read input VCF files ----------------
  
  vcf_list = lapply(vcf_list, function(x) {
    
    out = cbind(x@fix[,1:7],
                vcfR::extract_gt_tidy(x, verbose = FALSE),
                vcfR::extract_info_tidy(x))
    
    out = data.table::as.data.table(out)
    
    ANN_matrix = str_split(out$ANN, "\\|", simplify = TRUE)[,c(4, 10, 11)]
    # ANN_matrix = paste0(ANN_matrix[,1], "|", ANN_matrix[,2])
    
    out = out[,c("CHROM", 
                 "POS",
                 "ID",
                 "REF",
                 "ALT",
                 "gt_DP",
                 "gt_AD",
                 "TYPE"), with = FALSE]
    
    out$Gene_Name = ANN_matrix[,1]
    out$HGVS.c = ANN_matrix[,2]
    out$HGVS.p = ANN_matrix[,3]
    
    return(out)
    
  })
  
  # Covert VCF to table ------------------
  
  vcf_list = lapply(vcf_list, function(x) {
    
    ALT_matrix = str_split(x$ALT, ",", simplify = TRUE)
    gt_AD_matrix = str_split(x$gt_AD, ",", simplify = TRUE)
    TYPE_matrix = str_split(x$TYPE, ",", simplify = TRUE)
    
    out = list()
    
    for(i in 1:ncol(ALT_matrix)) {
      
      who = which(ALT_matrix[,i] != "")
      
      
      out[[i]] = data.table(CHROM = x[who, ]$CHROM,
                            POS = x[who, ]$POS,
                            ID = x[who, ]$ID,
                            REF = x[who, ]$REF,
                            ALT = ALT_matrix[who, i],
                            gt_DP = x[who, ]$gt_DP,
                            gt_AD_ref = gt_AD_matrix[who, 1],
                            gt_AD_alt = gt_AD_matrix[who, i+1],
                            TYPE = TYPE_matrix[who, i],
                            Gene_Name = x[who, ]$Gene_Name,
                            HGVS.c = x[who, ]$HGVS.c,
                            HGVS.p = x[who, ]$HGVS.p)
      
    }
    
    out = rbindlist(out)
    
    out$POS = as.numeric(out$POS)
    
    out = out[order(POS), ]
    
    return(out)
    
  })
  
  vcf_list = lapply(vcf_list, function(x) {
    
    x$gt_AD_ref = as.numeric(x$gt_AD_ref)
    x$gt_AD_alt = as.numeric(x$gt_AD_alt)
    
    x$AF = x$gt_AD_alt / x$gt_DP
    
    x$ID = paste(x$CHROM, x$POS, x$REF, x$ALT, sep = ";")
    
    return(x)
    
  })
  
  for(i in names(vcf_list)) {
    
    vcf_list[[i]]$sample = i
    
  }
  
  # Combine all vcfs into one table --------------
  
  vcf_list = rbindlist(vcf_list)
  
  vcf_list_collapsed = vcf_list[, .(n = .N,
                                    samples_found = paste(sample, collapse = "+")), by = ID]
  
  who = match(vcf_list$ID, vcf_list_collapsed$ID)
  vcf_list$samples_occur = vcf_list_collapsed[who, ]$samples_found
  
  vcf_list$freq = as.character(as.numeric(f) / 100)
  
  total_vcf_list[[f]] = vcf_list
  
}

total_vcf_list = rbindlist(total_vcf_list)

gr = ggplot(data = total_vcf_list, aes(x = samples_occur)) +
  
  geom_bar(aes(fill = samples_occur), color = "black", alpha = 0.7) +
  
  scale_fill_nejm() +
  
  theme_minimal() +
  
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 25, hjust = 1),
        panel.background = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA)) +
  
  facet_wrap(vars(freq), nrow = 1, scales = "free")


saveImageHigh::save_as_pdf({print(gr)},
                           file.name = "Nvariants-per-minFreq.pdf",
                           width = 12,
                           height = 3)


plot_matrix = total_vcf_list[, c("ID", "AF", "samples_occur", "sample", "freq"), with = FALSE]
plot_matrix = plot_matrix %>% tidyr::gather(key = "type", value = "value", -ID, -samples_occur, -sample, -freq)

B = ggplot(data = plot_matrix, aes(x = samples_occur, y = value, fill = sample)) +
  geom_boxplot(width = 0.3, alpha = 0.7) +
  scale_fill_nejm() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none") +
  facet_wrap(vars(freq), nrow = 1, scales = "free") 

saveImageHigh::save_as_pdf({print(B)},
                           file.name = "AF-per-minFreq.pdf",
                           width = 12,
                           height = 3)


vcf_list = total_vcf_list[which(total_vcf_list$freq == 0.1), ]

plot_matrix = vcf_list[, c("ID", "AF", "gt_DP", "gt_AD_alt", "sample"), with = FALSE]
plot_matrix = plot_matrix %>% tidyr::gather(key = "type", value = "value", -ID, -sample)

table(vcf_list$sample)

# Check mean 
t.test(x = plot_matrix[which(plot_matrix$sample == "Replicate_A" & plot_matrix$type == "gt_DP"), ]$value, 
       y = plot_matrix[which(plot_matrix$sample == "Replicate_B" & plot_matrix$type == "gt_DP"), ]$value)

t.test(x = plot_matrix[which(plot_matrix$sample == "Replicate_A" & plot_matrix$type == "gt_AD_alt"), ]$value, 
       y = plot_matrix[which(plot_matrix$sample == "Replicate_B" & plot_matrix$type == "gt_AD_alt"), ]$value)

t.test(x = plot_matrix[which(plot_matrix$sample == "Replicate_A" & plot_matrix$type == "AF"), ]$value, 
       y = plot_matrix[which(plot_matrix$sample == "Replicate_B" & plot_matrix$type == "AF"), ]$value)


# Fisher F test
var.test(x = plot_matrix[which(plot_matrix$sample == "Replicate_A" & plot_matrix$type == "gt_DP"), ]$value, 
         y = plot_matrix[which(plot_matrix$sample == "Replicate_B" & plot_matrix$type == "gt_DP"), ]$value)

var.test(x = plot_matrix[which(plot_matrix$sample == "Replicate_A" & plot_matrix$type == "gt_AD_alt"), ]$value, 
         y = plot_matrix[which(plot_matrix$sample == "Replicate_B" & plot_matrix$type == "gt_AD_alt"), ]$value)

var.test(x = plot_matrix[which(plot_matrix$sample == "Replicate_A" & plot_matrix$type == "AF"), ]$value, 
         y = plot_matrix[which(plot_matrix$sample == "Replicate_B" & plot_matrix$type == "AF"), ]$value)

