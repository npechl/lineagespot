

rm(list = ls())

# load libraries ---------------------------------------------------------------

source("libraries.R")
source("helpers.R")

# merge VCFs -------------------------------------------------------------------

vcf_list = merge_VCFs(vcfs_path = "../Thessaloniki/", 
                      metaData_path = "../Thessaloniki/meta-data.xlsx")

write.csv(vcf_list, 
          file = "sewage_VCFList.csv", 
          row.names = FALSE, 
          quote = FALSE)








# # Original DP, AD, AF statistics -----------------------------------------------
# 
# plot_matrix = vcf_list[, c("ID", "gt_DP", "gt_AD_alt", "AF", "sample"), with = FALSE]
# 
# ## number of reads - Pearson correlation  --------------------------------------
# 
# plot_matrix$nreads = 0
# 
# who = base::match(plot_matrix$sample, meta_data$Sample_ID)
# 
# plot_matrix$nreads = meta_data[who, ]$input_freebayes
# 
# cor.test(x = plot_matrix$gt_DP,
#          y = plot_matrix$nreads,
#          method = c("pearson"))
# 
# cor.test(x = plot_matrix$gt_AD_alt,
#          y = plot_matrix$nreads,
#          method = c("pearson"))
# 
# cor.test(x = plot_matrix$AF,
#          y = plot_matrix$nreads,
#          method = c("pearson"))

# number of variants per sample -----------------------------------------------

table(vcf_list$sample)


# mean and variance per sample ------------------------------------------------

plot_matrix = vcf_list[, c("ID", "gt_DP", "gt_AD_alt", "AF", "sample"), with = FALSE]

plot_matrix = plot_matrix %>% tidyr::gather(key = "type", value = "value", -ID, -sample)


raw_stats_matrix = compute_statistics(data = plot_matrix, samplelist = unique(plot_matrix$sample))


gr = ggplot(data = plot_matrix, aes(x = sample, y = value, fill = sample)) +
  
  geom_boxplot(width = 0.3, alpha = 0.7) +
  
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "none") +
  
  facet_wrap(vars(type), 
             ncol = 3, 
             nrow = 1, 
             scales = "free")


saveImageHigh::save_image({print(gr)},
                          file.name = "raw-stats.jpeg",
                          width = 8,
                          height = 6)



# check common variants -----------------------

vcf_list = common_variants(data = vcf_list)

common_stats = vcf_list$stats

vcf_list = vcf_list$data


plot_matrix = vcf_list[, c("ID", "gt_DP", "gt_AD_alt", "AF", "samples_occur", "sample"), with = FALSE]
plot_matrix = plot_matrix %>% tidyr::gather(key = "type", value = "value", -ID, -samples_occur, -sample)

gr = ggplot(data = plot_matrix, aes(x = samples_occur, y = value, fill = sample)) +
  geom_boxplot(width = 0.3, alpha = 0.7) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none") +
  facet_wrap(vars(type), ncol = 3, nrow = 1, scales = "free")

saveImageHigh::save_image({print(gr)},
                          file.name = "raw-common-stats_v2.jpeg",
                          width = 18,
                          height = 9)

