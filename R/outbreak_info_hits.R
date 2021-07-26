#' outbreak_info_hits
#'
#' Find overlapping variants with SARS-CoV-2 reference lineages
#' coming from outbreak.info reports
#'
#' @param vcf_table
#'
#' A tab-delimited table containing all variants for all samples
#'
#' @param ref_folder
#'
#' A path to outbreak.info reports
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
#' @export outbreak_info_hits
#'
#' @examples


outbreak_info_hits <- function(vcf_table = NULL,
                               ref_folder = NULL,
                               print.out = FALSE) {


  if( base::is.null(vcf_table) ) {

    base::stop('Please provide a tab-delimited table containing variants.')

  }

  if( base::is.null(ref_folder) ) {

    base::stop('Please provide a path to outbreak.info VoC reports.')

  }

  # create AA abbreviation --------------------------------------------------------

  AA_abbreviations = data.table::data.table(Three_Letter = c("Ala",
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

  vcf_table$AA_alt = stringr::str_remove(vcf_table$AA_alt, "p\\.")


  # clean AA variants ------------------------------------------------------------

  aa_variants = stringr::str_split(vcf_table$AA_alt, "[0-9]|_")
  codon_num = stringr::str_split(vcf_table$AA_alt, "[A-Z]|[a-z]|\\_|\\*|\\?")

  aa_variants = base::lapply(aa_variants, function(x) {

    return( x[x != ""] )

  })

  codon_num = base::lapply(codon_num, function(x) {

    return( x[x != ""] )

  })


  # merge overall tables ---------------------------------------------------------

  aa_split_list = base::lapply(seq(1:base::nrow(vcf_table)), function(i) {

    out = data.table::data.table(ref_aa = NA,
                                 alt_aa = NA,
                                 info = NA,
                                 codon_num = NA,
                                 codon_end = NA)


    if(base::length(aa_variants[[i]]) == 2) {

      out$ref_aa = aa_variants[[i]][1]
      out$alt_aa = aa_variants[[i]][2]

      out$codon_num = codon_num[[i]][1]

    } else if(base::length(aa_variants[[i]] == 3)) {

      out$ref_aa = aa_variants[[i]][1]
      out$alt_aa = aa_variants[[i]][2]

      out$type == aa_variants[[i]][3]

      out$codon_num = codon_num[[i]][1]
      out$codon_end = codon_num[[i]][2]

    } else {


    }

    return(out)

  })

  aa_split_list = data.table::rbindlist(aa_split_list)

  vcf_table = base::cbind(vcf_table, aa_split_list)

  vcf_table$change_length_nt = base::abs(stringr::str_length(vcf_table$ALT) - stringr::str_length(vcf_table$REF))

  vcf_table$codon_num = base::as.numeric(vcf_table$codon_num)
  vcf_table$codon_end = base::as.numeric(vcf_table$codon_end)


  # add AA variants with Nt when there is not any --------------------------------

  vcf_table[base::which(vcf_table$AA_alt == ""), ]$AA_alt = vcf_table[base::which(vcf_table$AA_alt == ""), ]$ID


  # clean AA abbreviation --------------------------------------------------------

  for(i in 1:base::nrow(AA_abbreviations)) {

    vcf_table$ref_aa = stringr::str_replace_all(vcf_table$ref_aa,
                                                AA_abbreviations[i,]$Three_Letter,
                                                AA_abbreviations[i,]$One_Letter)

    vcf_table$alt_aa = stringr::str_replace_all(vcf_table$alt_aa,
                                                AA_abbreviations[i,]$Three_Letter,
                                                AA_abbreviations[i,]$One_Letter)

    vcf_table$AA_alt = stringr::str_replace_all(vcf_table$AA_alt,
                                                AA_abbreviations[i,]$Three_Letter,
                                                AA_abbreviations[i,]$One_Letter)

  }


  # read reference ---------------------------------------------------------------

  fls_ref = base::list.files(ref_folder,
                             pattern = "outbreakinfo",
                             full.names = TRUE)

  VoC_hits_list = base::list()

  for(ref_index in fls_ref) {

    ref = data.table::fread(ref_index)

    ref[base::which(ref$gene == "ORF1a"), ] $gene = "ORF1ab"
    ref[base::which(ref$gene == "ORF1b"), ] $gene = "ORF1ab"

    ref = ref[base::order(ref$gene), ]


    ref$barcode = base::paste0(ref$gene, ":", ref$ref_aa, ref$codon_num, ref$alt_aa)

    ref[base::which(ref$type == "deletion"), ]$barcode = ref[base::which(ref$type == "deletion"), ]$ref_aa

    ref[base::which(ref$change_length_nt == "None"), ]$change_length_nt = "0"

    ref$change_length_nt = base::as.numeric(ref$change_length_nt)


    # VoC hits -----------------------------------------------------------------

    ref$alt_aa = stringr::str_replace(ref$alt_aa, "\\*", "\\\\*")
    vcf_table$alt_aa = stringr::str_replace(vcf_table$alt_aa, "\\*", "\\\\*")

    voc_data = base::list()

    for(i in 1:base::nrow(ref)) {


      if(ref[i,]$type == "substitution") {

        temp = vcf_table[base::which(vcf_table$Gene_Name == ref[i,]$gene), ]

        temp = temp[base::which(stringr::str_detect(temp$ref_aa, ref[i,]$ref_aa)), ]

        temp = temp[base::which(stringr::str_detect(temp$alt_aa, ref[i,]$alt_aa)), ]

        temp = temp[base::which(temp$codon_num == ref[i,]$codon_num | temp$codon_num == ref[i,]$codon_num + 1 | temp$codon_num == ref[i,]$codon_num - 1), ]

        temp = temp[base::which(temp$change_length_nt == ref[i,]$change_length_nt), ]

      } else {

        temp = vcf_table[base::which(vcf_table$Gene_Name == ref[i,]$gene), ]
        temp = temp[base::which(stringr::str_detect(temp$AA_alt, "del")), ]

        temp = temp[base::which(temp$codon_num == ref[i,]$codon_num), ]
        temp = temp[base::which(temp$codon_end == ref[i,]$codon_end), ]

        temp = temp[base::which(temp$change_length_nt == ref[i,]$change_length_nt), ]


      }


      voc_data[[ref[i,]$barcode]] = temp



    }

    not_overlapping_variants = base::names(which(lapply(voc_data, nrow) == 0))


    voc_data = data.table::rbindlist(voc_data)

    voc_data = base::unique(voc_data)



    # collapse table ---------------------------------------------------------------

    voc_data = voc_data[, c("POS",
                            "DP",
                            "AD_alt",
                            "Gene_Name",
                            "AA_alt",
                            "sample"), with = FALSE]

    voc_data = voc_data[, .(DP = unique(DP), AD_alt = sum(AD_alt)), by = .(POS, Gene_Name, AA_alt, sample)]

    # voc_data = voc_data[, base::list("DP" = base::unique(.SD[1]),
    #                                  "AD_alt" = base::sum(.SD[2])),
    #
    #                     by = base::list("POS", "Gene_Name", "AA_alt", "sample"),
    #
    #                     .SDcols = c("DP", "AD_alt")]

    # voc_positions = base::unique(voc_data[, c("POS", "AA_alt", "Gene_Name"), with = FALSE])

    # voc_data = voc_data[, base::list("DP" = base::sum(.SD[1]),
    #                                  "AD_alt" = base::sum(.SD[2])),
    #
    #                     by = base::list("Gene_Name", "AA_alt", "sample"),
    #
    #                     .SDcols = c("DP", "AD_alt")]

    voc_data = voc_data[, .(DP = sum(DP), AD_alt = sum(AD_alt)), by = .(Gene_Name, AA_alt, sample)]


    voc_data$AF = voc_data$AD_alt / voc_data$DP

    # add non overlapping rules ----------------------------------------------------


    not_overlapping_variants = str_split(not_overlapping_variants, "\\:", simplify = TRUE)

    voc_data = rbind(voc_data,
                     data.table(Gene_Name = rep(not_overlapping_variants[,1], length(unique(voc_data$sample))),
                                AA_alt = rep(not_overlapping_variants[,2], length(unique(voc_data$sample))),
                                sample = unique(voc_data$sample),
                                DP = 0,
                                AD_alt = 0,
                                AF = 0))

    VoC_hits_list[[ref_index]] = voc_data

  }


  # rm(list = setdiff(ls(), c("VoC_hits_list", "meta")))

  for(i in base::names(VoC_hits_list)) {

    strain = stringr::str_split(i, "_", simplify = TRUE)
    strain = strain[, base::ncol(strain)]
    strain = stringr::str_remove(strain, "\\.tsv")

    VoC_hits_list[[i]]$lineage = strain

  }

  # rm(i, strain)


  voc_data = data.table::rbindlist(VoC_hits_list)


  if( print.out ) {

    data.table::fwrite(voc_data,
                       file = "sewage_outbreak.info_hits.txt",
                       row.names = FALSE,
                       quote = FALSE,
                       sep = "\t")

  }

  base::return(voc_data)

}
