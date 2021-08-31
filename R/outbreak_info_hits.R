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
#' @param file.out
#'
#' Given name for the output file
#'
#' @import data.table
#' @import stringr
#'
#' @return
#'
#' @export lineagespot_hits
#'
#' @examples


lineagespot_hits <- function(vcf_table = NULL,
                             ref_folder = NULL,
                             file.out = "variantsList.lineagespot.txt",
                             print.out = FALSE) {


  if( is.null(vcf_table) ) {

    stop('Please provide a tab-delimited table containing variants.')

  }

  if( is.null(ref_folder) ) {

    stop('Please provide a path to outbreak.info VoC reports.')

  }

  # clean AA variants ------------------------------------------------------------

  aa_variants = str_split(vcf_table$AA_alt, "[0-9]|_")
  codon_num = str_split(vcf_table$AA_alt, "[A-Z]|[a-z]|\\_|\\*|\\?")

  aa_variants = base::lapply(aa_variants, function(x) {

    return( x[x != ""] )

  })

  codon_num = base::lapply(codon_num, function(x) {

    return( x[x != ""] )

  })


  # merge overall tables ---------------------------------------------------------

  aa_split_list = base::lapply(seq(1:nrow(vcf_table)), function(i) {

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

  vcf_table = cbind(vcf_table, aa_split_list)

  vcf_table$change_length_nt = base::abs(str_length(vcf_table$ALT) - str_length(vcf_table$REF))

  vcf_table$codon_num = as.numeric(vcf_table$codon_num)
  vcf_table$codon_end = as.numeric(vcf_table$codon_end)


  # add AA variants with Nt when there is not any --------------------------------

  # vcf_table[base::which(vcf_table$AA_alt == ""), ]$AA_alt = vcf_table[base::which(vcf_table$AA_alt == ""), ]$ID

  # read reference ---------------------------------------------------------------

  fls_ref = base::list.files(ref_folder,
                             pattern = "txt",
                             full.names = TRUE)

  VoC_hits_list = base::list()

  for(ref_index in fls_ref) {

    ref = fread(ref_index)

    ref = ref[order(ref$gene), ]

    ref$barcode = paste0(ref$gene, ":", ref$`amino acid`)

    ref$type = "snp"

    ref[which(str_detect(ref$`amino acid`, "del")), ]$type = "deletion"

    aa_variants = str_split(ref$`amino acid`, "[0-9]|_|\\/")
    codon_num = str_split(ref$`amino acid`, "[A-Z]|[a-z]|\\_|\\*|\\?|\\/")

    aa_variants = base::lapply(aa_variants, function(x) {

      x = x[x != ""]

      x = str_replace_all(x, "\\*", "\\\\*")

      return( x )

    })

    codon_num = base::lapply(codon_num, function(x) {

      return( as.numeric(x[x != ""]) )

    })

    # VoC hits -----------------------------------------------------------------

    # ref$alt_aa = stringr::str_replace(ref$alt_aa, "\\*", "\\\\*")
    # vcf_table$alt_aa = stringr::str_replace(vcf_table$alt_aa, "\\*", "\\\\*")

    voc_data = list()

    for(i in 1:nrow(ref)) {

      if(ref[i,]$type == "snp") {

        temp = vcf_table[which(vcf_table$Gene_Name == ref[i,]$gene), ]

        temp = temp[which(str_detect(temp$ref_aa, aa_variants[[i]][1] )), ]

        temp = temp[base::which(stringr::str_detect(temp$alt_aa, aa_variants[[i]][2] )), ]

        who = which(temp$codon_num == codon_num[[i]] | temp$codon_num == codon_num[[i]] + 1 | temp$codon_num == codon_num[[i]] - 1)

        temp = temp[who, ]

        # temp = temp[base::which(temp$change_length_nt == ref[i,]$change_length_nt), ]

      } else {

        temp = vcf_table[which(vcf_table$Gene_Name == ref[i,]$gene), ]
        temp = temp[which(str_detect(temp$AA_alt, "del")), ]

        who = which(temp$codon_num == codon_num[[i]][1] | temp$codon_num == codon_num[[i]][2] | temp$codon_end == codon_num[[i]][2] | temp$codon_end == codon_num[[i]][1])
        temp = temp[who, ]

        # temp = temp[which(temp$codon_end == codon_num[[i]][2]), ]
        # temp = temp[base::which(temp$change_length_nt == ref[i,]$change_length_nt), ]

      }


      voc_data[[ref[i,]$barcode]] = temp



    }

    not_overlapping_variants = names(which(lapply(voc_data, nrow) == 0))


    voc_data = rbindlist(voc_data)

    voc_data = unique(voc_data)



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

    strain = str_split(i, "\\/", simplify = TRUE)
    strain = strain[, ncol(strain)]
    strain = str_remove(strain, "\\.txt")

    VoC_hits_list[[i]]$lineage = strain

  }

  # rm(i, strain)


  voc_data = rbindlist(VoC_hits_list)


  if( print.out ) {

    fwrite(voc_data,
           file = file.out,
           row.names = FALSE,
           quote = FALSE,
           sep = "\t")

  }

  return(voc_data)

}
