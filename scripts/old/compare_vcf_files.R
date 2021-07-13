library(vcfR)
library(xlsx)

rm(list = ls())

vcf_list = list()

vcf_list[["s100"]] = read.vcfR("vcf-files/Sewage_CoV19_L3_S1_s100_freebayes.vcf", verbose = FALSE)
vcf_list[["s200"]] = read.vcfR("vcf-files/Sewage_CoV19_L3_S1_s200_freebayes.vcf", verbose = FALSE)
vcf_list[["s500"]] = read.vcfR("vcf-files/Sewage_CoV19_L3_S1_s500_freebayes.vcf", verbose = FALSE)

out = matrix(data = "0", nrow = 3, ncol = 3)

colnames(out) = c("Sewage_CoV19_L3_S1_s100", "Sewage_CoV19_L3_S1_s200", "Sewage_CoV19_L3_S1_s500")
row.names(out) = c("Sewage_CoV19_L3_S1_s100", "Sewage_CoV19_L3_S1_s200", "Sewage_CoV19_L3_S1_s500")

AF_comparison = list()

for(i in 1:3) {
  
  temp_i = paste0(vcf_list[[i]]@fix[,"POS"], ";", vcf_list[[i]]@fix[,"ALT"])
  
  AO = vcfR::extract_info_tidy(vcf_list[[i]])$AO
  AO = str_split(AO, ",", simplify = TRUE)
  
  AO = apply(AO, 2, as.numeric)
  
  DP = vcfR::extract_info_tidy(vcf_list[[i]])$DP
  
  AF_i = AO / DP
  
  AF_i = round(AF_i, digits = 2)
  
  for(j in 1:3) {
    
    temp_j = paste0(vcf_list[[j]]@fix[,"POS"], ";", vcf_list[[j]]@fix[,"ALT"])
    
    AO = vcfR::extract_info_tidy(vcf_list[[j]])$AO
    AO = str_split(AO, ",", simplify = TRUE)
    
    AO = apply(AO, 2, as.numeric)
    
    DP = vcfR::extract_info_tidy(vcf_list[[j]])$DP
    
    AF_j = AO / DP
    
    AF_j = round(AF_j, digits = 2)
    
    out[i, j] = paste0(length(which(temp_i %in% temp_j)), " ; ", 
                       min(AF_i[which(temp_i %in% temp_j), ], na.rm = TRUE), " ; ", 
                       round(mean(AF_i[which(temp_i %in% temp_j), ], na.rm = TRUE), digits = 2), " ; ", 
                       max(AF_i[which(temp_i %in% temp_j), ], na.rm = TRUE), 
                       " / ", 
                       length(which(!(temp_i %in% temp_j))), " ; ", 
                       min(AF_i[which(!(temp_i %in% temp_j)), ], na.rm = TRUE), " ; ", 
                       round(mean(AF_i[which(!(temp_i %in% temp_j)), ], na.rm = TRUE), digits = 2), " ; ", 
                       max(AF_i[which(!(temp_i %in% temp_j)), ], na.rm = TRUE))
    
    if(i < j) {
      
      AF_i_temp = AF_i[which(temp_i %in% temp_j), ]
      AF_j_temp = AF_j[which(temp_j %in% temp_i), ]
      
      # AF_i[is.na(AF_i)] = ""
      # AF_j[is.na(AF_j)] = ""
      
      AF_i_temp = apply(AF_i_temp, 1 , paste, collapse = ",")
      AF_j_temp = apply(AF_j_temp, 1 , paste, collapse = ",")
      
      AF_i_temp = str_remove_all(AF_i_temp, ",NA")
      AF_j_temp = str_remove_all(AF_j_temp, ",NA")
      
      AF_comparison[[paste0(i, "-", j)]] = data.frame(Mut = temp_i[which(temp_i %in% temp_j)],
                                                      AF1 = AF_i_temp,
                                                      AF2 = AF_j_temp)
      
    }
    
  }
  
}

write.xlsx(out, file = "size-dataset-comparison.xlsx", sheetName = "mutation-comp")

for(i in 1:length(AF_comparison)) {
  
  write.xlsx(AF_comparison[[i]], 
             file = "size-dataset-comparison.xlsx", 
             row.names = FALSE,
             sheetName = names(AF_comparison[i]),
             append = TRUE)
  
}




