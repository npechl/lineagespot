#' lineagespot function
#'
#'
#' Identify SARS-CoV-2 related mutations based on a single (or a list) of 
#' variant(s) file(s)
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
#' @param file.out.index
#'
#' A string index merged at the output files
#'
#' @param gff3_path
#'
#' Path to gff file
#'
#' @param ref_folder
#'
#' A path to outbreak.info reports
#'
#' @param voc
#'
#' a character vector containing the names of the lineages of interest
#'
#' @param AF_threshold
#'
#' A parameter indicating the AF threshold for identifying variants per sample
#'
#' @import jsonlite
#' @import httr
#' @import data.table
#' @import stringr
#' @import vcfR
#'
#' @return List
#'
#' @export merge_vcf
#'
#' @examples


lineagespot <- function(vcf_fls = NULL,
                        vcf_folder = NULL,
                        gff3_path = NULL,
                        ref_folder = NULL,
                        voc = c("B.1.617.2", "B.1.1.7", "B.1.351", "P.1"),
                        AF_threshold = 0.8,
                        file.out.index = Sys.Date(),
                        print.out = FALSE) {


  vcf_table = merge_vcf(vcf_fls = vcf_fls,
                        vcf_folder = vcf_folder,
                        gff3_path = gff3_path,
                        file.out = paste0("Variant_table_", file.out.index, ".txt"),
                        print.out = print.out)



  hits_table = lineagespot_hits(vcf_table = vcf_table,
                                ref_folder = ref_folder,
                                voc = voc,
                                file.out = paste0("lineage_hits_", file.out.index, ".txt"),
                                print.out = print.out)



  lineage_report = uniq_variants(hits_table = hits_table,
                                 AF_threshold = AF_threshold,
                                 file.out = paste0("lineage_report_", file.out.index, ".txt"),
                                 print.out = print.out)

  out = list("variants.table" = vcf_table,
             "lineage.hits" = hits_table,
             "lineage.report" = lineage_report)


  return(out)

}

# Removes R CMD check NOTE regarding global variables
utils::globalVariables(c("DP", ".", "AD_alt", "POS", "Gene_Name",
    "AA_alt", "start_pos", "end_pos", "gene_name", "lineage", 
    "AF", ""))











