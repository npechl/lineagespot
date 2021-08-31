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
#' @param file.out
#'
#' Given name for the output file
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
                      file.out = "vcfList_table.txt",
                      print.out = FALSE) {


  if( is.null(vcf_fls) & is.null(vcf_folder) ) {

    stop("Please provide some VCF files")

  }


  if( is.null(vcf_fls) ) {

    vcf_fls = list.files(vcf_folder,
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

  # add variant parameters ----------------------------------------------------------

  vcf_list = base::lapply(vcf_list, change_AA_abbreviations)

  # add variant parameters ----------------------------------------------------------

  vcf_list = base::lapply(vcf_list, correct_Orf1ab_gene)

  # add sample name ------------------------------------------------------------------

  for(i in names(vcf_list)) {

    sample_name = base::unlist(str_split(i, "\\/"))
    sample_name = sample_name[length(sample_name)]

    sample_name = str_split(sample_name, "\\.vcf", simplify = TRUE)[,1]

    vcf_list[[i]]$sample = sample_name

  }

  # rm(i)
  # gc()

  # Combine all vcfs into one table -----------------------------------------------------

  vcf_list = rbindlist(vcf_list)


  if( print.out ) {

    fwrite(vcf_list,
           file = file.out,
           row.names = FALSE,
           quote = FALSE,
           sep = "\t")

  }

  return(vcf_list)

}

vcf_to_table <- function(x) {

  out = cbind(x@fix[,1:7],
              vcfR::extract_gt_tidy(x, verbose = FALSE),
              vcfR::extract_info_tidy(x))

  out = as.data.table(out)

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

  out = out[order(out$POS), ]

  return(out)

}

compute_AF <- function(x) {

  x$AD_ref = as.numeric(x$AD_ref)
  x$AD_alt = as.numeric(x$AD_alt)

  x$AF = x$AD_alt / x$DP

  x$ID = paste(x$CHROM, x$POS, x$REF, x$ALT, sep = ";")

  return(x)

}

change_AA_abbreviations <- function(x) {

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


  x$Nt_alt = str_remove_all(x$Nt_alt, "c\\.")
  x$AA_alt = str_remove_all(x$AA_alt, "p\\.")

  for(i in 1:nrow(AA_abbreviations)) {

    x$AA_alt = str_replace_all(x$AA_alt,
                               AA_abbreviations[i,]$Three_Letter,
                               AA_abbreviations[i,]$One_Letter)

  }

  return(x)

}

correct_Orf1ab_gene <- function(x) {

  # print(head(x))

  genes = data.table(gene_name = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"),

                     start_pos = c(266, 13442, 21563, 25393, 26245, 26523, 27202, 27394, 27756, 27894, 28274, 29558),

                     end_pos = c(13441, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 29533, 29674))


  x$codon_num = 0

  x = x[order(x$POS), ]

  for(i in 1:nrow(genes)) {

    # x = x[which(x$POS >= genes[i,]$start_pos & x$POS <= genes[i,]$end_pos), ]

    # print(nrow(x))

    x[which(x$POS >= genes[i,]$start_pos & x$POS <= genes[i,]$end_pos), ]$codon_num = x[which(x$POS >= genes[i,]$start_pos & x$POS <= genes[i,]$end_pos), ]$POS - genes[i,]$start_pos + 1
    x[which(x$POS >= genes[i,]$start_pos & x$POS <= genes[i,]$end_pos), ]$Gene_Name = genes[i,]$gene_name
    x[which(x$POS >= genes[i,]$start_pos & x$POS <= genes[i,]$end_pos), ]$codon_num = ( x[which(x$POS >= genes[i,]$start_pos & x$POS <= genes[i,]$end_pos), ]$codon_num %/% 3 ) + ( x[which(x$POS >= genes[i,]$start_pos & x$POS <= genes[i,]$end_pos), ]$codon_num %% 3 > 0 )

  }

  x = x[which(x$codon_num != 0), ]

  return(x)

}

