
rm(list = ls())

gc()

# load libraries ---------------------------------------------------------------

library(data.table)
library(stringr)

# read input data --------------------------------------------------------------

# please provide a path to the folder containing VCF files
path_to_vcfList = "vcfList_table.txt"

path_to_outbreak_info = "ref/outbreak_info/"

sewage_data = data.table::fread(path_to_vcfList)


# create AA abbreviation --------------------------------------------------------

AA_abbreviations = data.table(Three_Letter = c("Ala",
                                               "Arg",
                                               "Asn",
                                               "Asp",
                                               "Cys",
                                               "Glu",
                                               "Gln",
                                               "Gly",
                                               "His",
                                               "Ile",
                                               "Leu",
                                               "Lys",
                                               "Met",
                                               "Phe",
                                               "Pro",
                                               "Ser",
                                               "Thr",
                                               "Trp",
                                               "Tyr",
                                               "Val"),
                              
                              One_Letter = c("A",
                                             "R",
                                             "N",
                                             "D",
                                             "C",
                                             "E",
                                             "Q",
                                             "G",
                                             "H",
                                             "I",
                                             "L",
                                             "K",
                                             "M",
                                             "F",
                                             "P",
                                             "S",
                                             "T",
                                             "W",
                                             "Y",
                                             "V"))

# remove "p." character and correct ORF ----------------------------------------

sewage_data$AA_alt = str_remove(sewage_data$AA_alt, "p\\.") 


# clean AA variants ------------------------------------------------------------

aa_variants = str_split(sewage_data$AA_alt, "[0-9]|_")
codon_num = str_split(sewage_data$AA_alt, "[A-Z]|[a-z]|\\_|\\*|\\?")

aa_variants = lapply(aa_variants, function(x) {
  
  return( x[x != ""] )
  
})

codon_num = lapply(codon_num, function(x) {
  
  return( x[x != ""] )
  
})


# merge overall tables ---------------------------------------------------------

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

sewage_data$change_length_nt = abs(str_length(sewage_data$ALT) - str_length(sewage_data$REF))

sewage_data$codon_num = as.numeric(sewage_data$codon_num)
sewage_data$codon_end = as.numeric(sewage_data$codon_end)


# add AA variants with Nt when there is not any --------------------------------

sewage_data[which(sewage_data$AA_alt == ""), ]$AA_alt = sewage_data[which(sewage_data$AA_alt == ""), ]$ID


# clean AA abbreviation --------------------------------------------------------

for(i in 1:nrow(AA_abbreviations)) {
  
  sewage_data$ref_aa = str_replace_all(sewage_data$ref_aa, 
                                       AA_abbreviations[i,]$Three_Letter,
                                       AA_abbreviations[i,]$One_Letter)
  
  sewage_data$alt_aa = str_replace_all(sewage_data$alt_aa, 
                                       AA_abbreviations[i,]$Three_Letter,
                                       AA_abbreviations[i,]$One_Letter)
  
  sewage_data$AA_alt = str_replace_all(sewage_data$AA_alt,
                                       AA_abbreviations[i,]$Three_Letter,
                                       AA_abbreviations[i,]$One_Letter)
  
}


# read reference ---------------------------------------------------------------

fls_ref = list.files(path_to_outbreak_info, 
                     pattern = "outbreakinfo", 
                     full.names = TRUE)

VoC_hits_list = list()

for(ref_index in fls_ref) {
  
  ref = fread(ref_index)
  
  ref[which(ref$gene == "ORF1a"), ] $gene = "ORF1ab"
  ref[which(ref$gene == "ORF1b"), ] $gene = "ORF1ab"
  
  ref = ref[order(ref$gene), ]
  
  
  ref$barcode = paste0(ref$gene, ":", ref$ref_aa, ref$codon_num, ref$alt_aa)
  
  ref[which(ref$type == "deletion"), ]$barcode = ref[which(ref$type == "deletion"), ]$ref_aa
  
  ref[which(ref$change_length_nt == "None"), ]$change_length_nt = "0"
  
  ref$change_length_nt = as.numeric(ref$change_length_nt)
  
  
  # VoC hits -----------------------------------------------------------------
  
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
      
      temp = temp[which(temp$codon_num == ref[i,]$codon_num | temp$codon_num == ref[i,]$codon_num + 1 | temp$codon_num == ref[i,]$codon_num - 1), ]
      
      temp = temp[which(temp$change_length_nt == ref[i,]$change_length_nt), ]
      
    } else {
      
      temp = sewage_data[which(sewage_data$Gene_Name == ref[i,]$gene), ]
      temp = temp[which(str_detect(temp$AA_alt, "del")), ]
      
      temp = temp[which(temp$codon_num == ref[i,]$codon_num), ]
      temp = temp[which(temp$codon_end == ref[i,]$codon_end), ]
      
      temp = temp[which(temp$change_length_nt == ref[i,]$change_length_nt), ]
      
      
    }
    
    
    voc_data[[ref[i,]$barcode]] = temp
    
    
    
  }
  
  not_overlapping_variants = names(which(lapply(voc_data, nrow) == 0))
  
  
  voc_data = rbindlist(voc_data)
  
  voc_data = unique(voc_data)
  
  
  
  # collapse table ---------------------------------------------------------------
  
  voc_data = voc_data[, c("POS", "DP", "AD_alt", "Gene_Name", "AA_alt", "sample"), 
                      with = FALSE]
  
  
  
  voc_data = voc_data[, .(DP = unique(DP),
                          AD_alt = sum(AD_alt)), 
                      by = .(POS, Gene_Name, AA_alt, sample)]
  
  voc_positions = unique(voc_data[, c("POS", "AA_alt", "Gene_Name"), with = FALSE])
  
  voc_data = voc_data[, .(DP = sum(DP),
                          AD_alt = sum(AD_alt)), 
                      by = .(Gene_Name, AA_alt, sample)]
  
  
  voc_data$AF = voc_data$AD_alt / voc_data$DP
  
  
  # add non overlapping rules ----------------------------------------------------
  
  
  # not_overlapping_variants = str_split(not_overlapping_variants, "\\:", simplify = TRUE)
  # 
  # voc_data = rbind(voc_data, 
  #                  data.table(Gene_Name = rep(not_overlapping_variants[,1], length(unique(voc_data$sample))),
  #                             HGVS.p = rep(not_overlapping_variants[,2], length(unique(voc_data$sample))),
  #                             sample = unique(voc_data$sample),
  #                             DP = 0,
  #                             AD = 0,
  #                             AF = 0))
  
  
  # rm(list = setdiff(ls(), c("voc_data", "meta")))
  
  
  VoC_hits_list[[ref_index]] = voc_data
  
}


rm(list = setdiff(ls(), c("VoC_hits_list", "meta")))

for(i in names(VoC_hits_list)) {
  
  strain = str_split(i, "_", simplify = TRUE)
  strain = strain[,ncol(strain)]
  strain = str_remove(strain, "\\.tsv")
  
  VoC_hits_list[[i]]$lineage = strain
  
}

rm(i, strain)


voc_data = rbindlist(VoC_hits_list)

data.table::fwrite(voc_data, 
                   file = "sewage_outbreak.info_hits.txt",
                   row.names = FALSE,
                   quote = FALSE,
                   sep = "\t")


# unique_lineage_variants = voc_data[, c("lineage", "HGVS.p"), with = FALSE]
# 
# unique_lineage_variants = unique(unique_lineage_variants)
# 
# unique_lineage_variants = names(which(table(unique_lineage_variants$HGVS.p) == 1))
# 
# voc_data = voc_data[which(voc_data$HGVS.p %in% unique_lineage_variants), ]


