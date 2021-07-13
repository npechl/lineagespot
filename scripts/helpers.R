# convert vcf format to table ------------------------------

vcf_to_table <- function(x) {
  
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

  colnames(out) = c("CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "DP",
                    "AD",
                    "TYPE")
  
  out$CHROM = str_squish(out$CHROM)
  out$Gene_Name = ANN_matrix[,1]
  out$Nt_alt = ANN_matrix[,2]
  out$AA_alt = ANN_matrix[,3]
  
  return(out)
  
}

# break multiple variants on same position -----------------

break_multiple_variants <- function(x) {
  
  ALT_matrix = str_split(x$ALT, ",", simplify = TRUE)
  AD_matrix = str_split(x$AD, ",", simplify = TRUE)
  TYPE_matrix = str_split(x$TYPE, ",", simplify = TRUE)
  
  out = list()
  
  for(i in 1:ncol(ALT_matrix)) {
    
    who = which(ALT_matrix[,i] != "")
    
    
    out[[i]] = data.table(CHROM = x[who, ]$CHROM,
                          POS = x[who, ]$POS,
                          ID = x[who, ]$ID,
                          REF = x[who, ]$REF,
                          ALT = ALT_matrix[who, i],
                          DP = x[who, ]$DP,
                          AD_ref = AD_matrix[who, 1],
                          AD_alt = AD_matrix[who, i+1],
                          TYPE = TYPE_matrix[who, i],
                          Gene_Name = x[who, ]$Gene_Name,
                          Nt_alt = x[who, ]$Nt_alt,
                          AA_alt = x[who, ]$AA_alt)
    
  }
  
  out = rbindlist(out)
  
  out$POS = as.numeric(out$POS)
  
  out = out[order(POS), ]
  
  return(out)
  
}

# add simple statistics ------------------------------------

compute_AF <- function(x) {
  
  x$AD_ref = as.numeric(x$AD_ref)
  x$AD_alt = as.numeric(x$AD_alt)
  
  x$AF = x$AD_alt / x$DP
  
  x$ID = paste(x$CHROM, x$POS, x$REF, x$ALT, sep = ";")
  
  return(x)
  
}
