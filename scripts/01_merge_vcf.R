
rm(list = ls())

# please provide a path to the folder containing VCF files
path_to_vcfFolder = "vcf-files/"

# please provide a tab-delimited meta data file
path_to_metaData = "vcf-files/meta-data.txt"

# load libraries ---------------------------------------------------------------

library(data.table)
library(stringr)

library(vcfR)

source("scripts/helpers.R")

# read meta data ---------------------------------------------------------------

meta_data = data.table::fread(path_to_metaData)

samplelist = meta_data[[1]]

# Read input VCF files ----------------------------------------------------------

vcf_list = list()

vcf_fls = base::list.files(path_to_vcfFolder, pattern = "vcf", full.names = TRUE)

for(i in samplelist) {
  
  vcf_list[[i]] = vcfR::read.vcfR(vcf_fls[which(str_detect(vcf_fls, i))], 
                                  verbose = FALSE)
  
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
  
  vcf_list[[i]]$sample = i
  
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