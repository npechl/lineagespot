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
#' @param gff3_path
#'
#' Path to gff file
#'
#' @param ref_folder
#'
#' A path to outbreak.info reports
#'
#' @param AF_threshold
#'
#' A parameter indicating the AF threshold for identifying variants per sample
#'
#' @import data.table
#' @import stringr
#'
#' @return
#'
#' @export merge_vcf
#'
#' @examples


lineagespot <- function(vcf_fls = NULL,
                        vcf_folder = NULL,
                        gff3_path = NULL,
                        ref_folder = NULL,
                        AF_threshold = 0.8,
                        file.out.index = "lineagespot",
                        print.out = FALSE) {


  vcf_table = merge_vcf(vcf_fls = vcf_fls,
                        vcf_folder = vcf_folder,
                        gff3_path = gff3_path,
                        file.out = paste0("vcfList_table_", file.out.index, ".txt"),
                        print.out = print.out)



  hits_table = lineagespot_hits(vcf_table = vcf_table,
                                ref_folder = ref_folder,
                                file.out = paste0("variantsList_", file.out.index, ".txt"),
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












