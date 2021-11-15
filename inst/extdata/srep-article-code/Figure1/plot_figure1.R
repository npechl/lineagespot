
rm(list = ls())

# load libraries ------------------------------------

library(data.table)
library(stringr)
library(ggplot2)

# Figure 1.A ---------------------------------------------


freebayes_gatk = fread("freebayes-vs-gatk.csv")
freebayes_gatk = freebayes_gatk[, -2]
freebayes_gatk$group = "freebayes vs gatk"

freebayes_mpileup = fread("freebayes-vs-mpileup.csv")
freebayes_mpileup = freebayes_mpileup[,-2]
freebayes_mpileup$group = "freebayes vs mpileup"

gatk_mpileup = fread("gatk-vs-mpileup.csv")
gatk_mpileup = gatk_mpileup[,-2]
gatk_mpileup$group = "gatk vs mpileup"


freebayes_gatk = freebayes_gatk[, 8:9]
freebayes_mpileup = freebayes_mpileup[, 8:9]
gatk_mpileup = gatk_mpileup[, 8:9]

colnames(freebayes_gatk) = c("diff", "group")
colnames(freebayes_mpileup) = c("diff", "group")
colnames(gatk_mpileup) = c("diff", "group")

plot_data = rbind(freebayes_gatk,
                  freebayes_mpileup,
                  gatk_mpileup)


gr1 = ggplot(data = plot_data, aes(x = diff, group = group)) +
  
  geom_density(adjust = 1,
               alpha = 0.75,
               fill = "gray50") +
  
  scale_y_continuous(limits = c(0, 2),
                     
                     expand = c(0, 0)) +
    
    scale_x_continuous(expand = c(0, 0)) +
  
  theme_minimal() +
  
  theme(# axis.title.x = element_blank(),
        
        # panel.border = element_rect(fill = NA),
        
        axis.ticks = element_line(),
        
        panel.grid = element_blank(),
        
        axis.line = element_line()
        ) +
  
  facet_wrap(~group, nrow = 1, scales = "free_y") + 
  
  labs(x = "Absolute difference of the number of common mutations between the three variant callers",
                                      y = "Density")


# Figure 1.B ------------------------------------------------

library(vcfR)

library(data.table)
library(stringr)

library(dplyr)

library(ggplot2)
library(ggsci)

# rm(list = ls())

# read meta data --------------------------------

meta_data = readxl::read_excel("input_files/meta-data.xlsx")

# Read sample list -----------------

# samplelist = read.table("paper_data/SampleList")
samplelist = meta_data[[1]]

freq_values = c("01", "10", "20", "30", "40", "50")

total_vcf_list = list()

for(f in freq_values) {
  
  vcf_list = list()
  
  for(i in samplelist) {
    
    vcf_list[[i]] = read.vcfR(paste0("input_files/", i, "_freebayes_", 
                                     f, "_ann.vcf"), verbose = FALSE)
    
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
  
  rm(i)
  gc()
  
  # Combine all vcfs into one table --------------
  
  vcf_list = rbindlist(vcf_list)
  
  vcf_list_collapsed = vcf_list %>% 
      group_by(ID) %>% 
      summarise(n = n(),  samples_found = paste(sample, collapse = "+"))
  
  
  who = match(vcf_list$ID, vcf_list_collapsed$ID)
  vcf_list$samples_occur = vcf_list_collapsed[who, ]$samples_found
  
  vcf_list$freq = as.character(as.numeric(f) / 100)
  
  total_vcf_list[[f]] = vcf_list
  
}

total_vcf_list = rbindlist(total_vcf_list)

total_vcf_list$freq = paste0("Low freq cutoff: ", total_vcf_list$freq)

gr2 = ggplot(data = total_vcf_list, aes(x = samples_occur)) + 
  geom_bar(aes(fill = samples_occur), color = "black", alpha = 0.5) +
  
  scale_fill_manual(values = c("gray50", "gray30", "gray70")) +
  
  # scale_fill_nejm() +
  
  # scale_fill_manual(values = c("Sample A" = "red4",
  #                              "Sample A+Sample B" = "cyan4",
  #                              "Sample B" = "blue4")) +
  
  theme_minimal()+
  
  scale_y_continuous(expand = c(0, 0 ), limits = c(0, 350)) +
  
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(),
        panel.background = element_blank(), 
        panel.grid = element_blank(),
        # panel.border = element_rect(fill = NA),
        axis.line = element_line()) +
  
  facet_wrap(vars(freq), nrow = 1, scales = "free") + 
  labs(y = "Number of Nt. mutations")


# Figure 1.C ----------------------------------------------------------

plot_matrix = total_vcf_list[, c("ID", "AF", "samples_occur", 
                                 "sample", "freq"), with = FALSE]

plot_matrix = plot_matrix %>% 
    tidyr::gather(key = "type", value = "value", 
                  -ID, -samples_occur, -sample, -freq)


gr3 = ggplot(data = plot_matrix, 
             aes(x = samples_occur, y = value, fill = sample)) +

  
  geom_point(aes(color = sample),
             alpha = 0.6,
             position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.3)) +
  
  geom_boxplot(width = 0.5, alpha = 0.4, outlier.shape = NA) +
  
  scale_color_manual(values = c("#00A087FF", "#F39B7FFF")) +
  scale_fill_manual(values = c("#00A087FF", "#F39B7FFF")) +
  
  theme_minimal() +
  theme(panel.grid = element_blank(),
        # panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(),
        axis.ticks = element_line(),
        legend.position = "none") +
  facet_wrap(vars(freq), nrow = 1, scales = "free") +
  
  ylim(0, 1) +
  
  labs(y = "Allele frequency")



library(patchwork)

gr1 / gr2 / gr3 + plot_annotation(tag_levels = "A")



# citation("data.table")
# citation("stringr")
# citation("ggplot2")
# citation("vcfR")
# citation("data.table")
# citation("stringr")
# citation("dplyr")
# citation("tidyr")
# citation("ggplot2")
# citation("ggsci")
# citation("patchwork")
