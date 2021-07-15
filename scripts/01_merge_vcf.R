
rm(list = ls())

# please provide a path to the folder containing VCF files
path_to_vcfFolder = "vcf-files/"

# load libraries ---------------------------------------------------------------

library(data.table)
library(stringr)

library(vcfR)

source("scripts/helpers.R")

# Read input VCF files ----------------------------------------------------------

vcf_list = list()

vcf_fls = base::list.files(path_to_vcfFolder, pattern = "vcf", full.names = TRUE)

for(i in vcf_fls) {
  
  vcf_list[[i]] = vcfR::read.vcfR(i, verbose = FALSE)
  
}

rm(i)


# Convert VCF to table ------------------------------------------------------------

vcf_list = lapply(vcf_list, vcf_to_table)

# break multiple rules ------------------------------------------------------------

vcf_list = lapply(vcf_list, break_multiple_variants)

# add variant parameters ----------------------------------------------------------

vcf_list = lapply(vcf_list, compute_AF)

# add sample name ------------------------------------------------------------------

for(i in names(vcf_list)) {
  
  sample_name = unlist(str_split(i, "\\/"))
  sample_name = sample_name[length(sample_name)]
  
  sample_name = str_split(sample_name, "_S", simplify = TRUE)[,1]
  
  vcf_list[[i]]$sample = sample_name
  
}

rm(i)
gc()

# Combine all vcfs into one table -----------------------------------------------------

vcf_list = rbindlist(vcf_list)

data.table::fwrite(vcf_list, 
                   file = "vcfList_table.txt",
                   row.names = FALSE,
                   quote = FALSE, 
                   sep = "\t")
