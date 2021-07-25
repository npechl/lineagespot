#' merge_vcf
#'
#' Merge Variant Calling Format (VCF) files into a single tab-delimited table
#'
#' @param vcf_fls
#'
#' A list of paths to vcf files
#'
#' @param vcf_folder
#'
#' A path to a folder containing all VCF file that will be integrated into a single table
#'
#' @param print.out
#'
#' logical value indicating if the produced table should be printed
#'
#' @import data.table
#' @import stringr
#'
#' @return
#'
#' @export merge_vcf
#'
#' @examples


merge_vcf <- function(vcf_fls = NULL,
                      vcf_folder = NULL,
                      print.out = FALSE) {


  if( base::is.null(vcf_fls) & base::is.null(vcf_folder) ) {

    base::stop("Please provide some VCF files")

  }


  if( base::is.null(vcf_fls) ) {

    vcf_fls = base::list.files(vcf_folder,
                               pattern = "vcf",
                               full.names = TRUE)


  }


  # Read input VCF files ----------------------------------------------------------

  vcf_list = base::list()

  for(i in vcf_fls) {

    vcf_list[[i]] = vcfR::read.vcfR(i, verbose = FALSE)

  }

  rm(i)


  # Convert VCF to table ------------------------------------------------------------

  vcf_list = base::lapply(vcf_list, vcf_to_table)

  # break multiple rules ------------------------------------------------------------

  vcf_list = base::lapply(vcf_list, break_multiple_variants)

  # add variant parameters ----------------------------------------------------------

  vcf_list = base::lapply(vcf_list, compute_AF)

  # add sample name ------------------------------------------------------------------

  for(i in names(vcf_list)) {

    sample_name = base::unlist(stringr::str_split(i, "\\/"))
    sample_name = sample_name[base::length(sample_name)]

    # sample_name = stringr::str_split(sample_name, "_S", simplify = TRUE)[,1]

    vcf_list[[i]]$sample = sample_name

  }

  # rm(i)
  # gc()

  # Combine all vcfs into one table -----------------------------------------------------

  vcf_list = data.table::rbindlist(vcf_list)


  if( print.out ) {

    data.table::fwrite(vcf_list,
                       file = "vcfList_table.txt",
                       row.names = FALSE,
                       quote = FALSE,
                       sep = "\t")

  }

  base::return(vcf_list)

}

vcf_to_table <- function(x) {

  out = base::cbind(x@fix[,1:7],
                    vcfR::extract_gt_tidy(x, verbose = FALSE),
                    vcfR::extract_info_tidy(x))

  out = data.table::as.data.table(out)

  ANN_matrix = stringr::str_split(out$ANN, "\\|", simplify = TRUE)[,c(4, 10, 11)]
  # ANN_matrix = paste0(ANN_matrix[,1], "|", ANN_matrix[,2])

  out = out[,c("CHROM",
               "POS",
               "ID",
               "REF",
               "ALT",
               "gt_DP",
               "gt_AD",
               "TYPE"), with = FALSE]

  base::colnames(out) = c("CHROM",
                          "POS",
                          "ID",
                          "REF",
                          "ALT",
                          "DP",
                          "AD",
                          "TYPE")

  out$CHROM = stringr::str_squish(out$CHROM)
  out$Gene_Name = ANN_matrix[,1]
  out$Nt_alt = ANN_matrix[,2]
  out$AA_alt = ANN_matrix[,3]

  return(out)

}

break_multiple_variants <- function(x) {

  ALT_matrix = stringr::str_split(x$ALT, ",", simplify = TRUE)
  AD_matrix = stringr::str_split(x$AD, ",", simplify = TRUE)
  TYPE_matrix = stringr::str_split(x$TYPE, ",", simplify = TRUE)

  out = list()

  for(i in 1:base::ncol(ALT_matrix)) {

    who = base::which(ALT_matrix[,i] != "")


    out[[i]] = data.table::data.table(CHROM = x[who, ]$CHROM,
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

  out = data.table::rbindlist(out)

  out$POS = base::as.numeric(out$POS)

  out = out[base::order(out$POS), ]

  return(out)

}

compute_AF <- function(x) {

  x$AD_ref = base::as.numeric(x$AD_ref)
  x$AD_alt = base::as.numeric(x$AD_alt)

  x$AF = x$AD_alt / x$DP

  x$ID = base::paste(x$CHROM, x$POS, x$REF, x$ALT, sep = ";")

  return(x)

}

