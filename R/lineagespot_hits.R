#' lineagespot_hits
#'
#' @description
#' Find overlapping variants with SARS-CoV-2 reference lineages
#' coming from outbreak.info reports
#'
#' @param vcf_table
#' A tab-delimited table containing all variants for all samples
#'
#' @param ref_folder
#' A path to lineages' reports
#'
#' @param voc
#' A character vector containing the names of the lineages of interest
#'
#' @param print.out
#' Logical value indicating if the produced table should be printed
#'
#' @param file.out
#' Given name for the output file
#'
#'
#' @return
#' A data table containing all identified SARS-CoV-2 variants
#' based on the provided reference files
#'
#' @export
#'
#' @examples
#'
#' variants_table = merge_vcf(vcf_folder = system.file("extdata",
#'                                                     "vcf-files",
#'                                                     package = "lineagespot"),
#'
#'                            gff3_path = system.file("extdata",
#'                                                    "NC_045512.2_annot.gff3",
#'                                                    package = "lineagespot"))
#'
#' # retrieve lineage reports using outbreak.info's API
#'
#' lineage_hits_table = lineagespot_hits(vcf_table = variants_table,
#'                                       voc = c("B.1.1.7", "B.1.617.2"))
#'
#'
#' # use user-specified references
#' lineage_hits_table.2 = lineagespot_hits(vcf_table = variants_table,
#'                                         ref_folder = system.file("extdata", "ref",
#'                                                                  package = "lineagespot"))
#'


lineagespot_hits <- function(vcf_table = NULL,
                             ref_folder = NULL,
                             voc = c("B.1.617.2", "B.1.1.7", "B.1.351", "P.1"),

                             file.out = paste0("lineage_hits_",
                                               Sys.Date(), ".txt"),

                             print.out = FALSE) {


    if( is.null(vcf_table) ) {

        stop('Please provide a tab-delimited table containing variants.')

    }

    # clean AA variants --------------------------------------------------------

    aa_variants = str_split(vcf_table$AA_alt, "[0-9]|_")
    codon_num = str_split(vcf_table$AA_alt, "[A-Z]|[a-z]|\\_|\\*|\\?")

    aa_variants = lapply(aa_variants, function(x) {

        return( x[x != ""] )

    })

    codon_num = lapply(codon_num, function(x) {

        return( x[x != ""] )

    })


    # merge overall tables -----------------------------------------------------

    aa_split_list = lapply(seq_len(nrow(vcf_table)), function(i) {

        out = data.table(ref_aa = NA,
                         alt_aa = NA,
                         info = NA,
                         codon_num = NA,
                         codon_start = NA,
                         codon_end = NA)


        if(length(aa_variants[[i]]) == 2) {

            out$ref_aa = aa_variants[[i]][1]
            out$alt_aa = aa_variants[[i]][2]

            out$codon_num = codon_num[[i]][1]

        } else if(length(aa_variants[[i]] == 3)) {

            out$ref_aa = aa_variants[[i]][1]
            out$alt_aa = aa_variants[[i]][2]

            out$type == aa_variants[[i]][3]

            out$codon_start = codon_num[[i]][1]
            out$codon_end = codon_num[[i]][2]

        } else {


        }

        return(out)

    })

    aa_split_list = rbindlist(aa_split_list)

    aa_split_list = aa_split_list[, c("ref_aa", "alt_aa", "info",
                                      "codon_start", "codon_end"), with = FALSE]

    vcf_table = cbind(vcf_table, aa_split_list)

    vcf_table$change_length_nt = abs(str_length(vcf_table$ALT) -
                                         str_length(vcf_table$REF))

    vcf_table$codon_num = as.numeric(vcf_table$codon_num)
    vcf_table$codon_start = as.numeric(vcf_table$codon_start)
    vcf_table$codon_end = as.numeric(vcf_table$codon_end)


    # read reference -----------------------------------------------------------

    if(!is.null(ref_folder)) {

        fls_ref = list.files(ref_folder,
        pattern = "txt",
        full.names = TRUE)

        reference_list = list()

        for(i in fls_ref){

            ref = fread(i)

            ref = ref[order(ref$gene), ]

            ref$barcode = paste0(ref$gene, ":", ref$`amino acid`)

            ref$type = "snp"

            ref[which(str_detect(ref$`amino acid`, "del")), ]$type = "deletion"

            reference_list[[i]] = ref

        }

    } else {

        reference_list = get_lineage_report(lineages = voc)

    }

    VoC_hits_list = list()

    for(ref_index in names(reference_list)) {

        ref = reference_list[[ref_index]]

        aa_variants = str_split(ref$`amino acid`, "[0-9]|_|\\/")
        codon_num = str_split(ref$`amino acid`, "[A-Z]|[a-z]|\\_|\\*|\\?|\\/")

        aa_variants = lapply(aa_variants, function(x) {

        x = x[x != ""]

        x = str_replace_all(x, "\\*", "\\\\*")

            return( x )

        })

        codon_num = lapply(codon_num, function(x) {

            return( as.numeric(x[x != ""]) )

        })

        # VoC hits -------------------------------------------------------------

        voc_data = list()

        for(i in seq_len(nrow(ref))) {

            if(ref[i,]$type == "snp") {

                temp = vcf_table[which(vcf_table$Gene_Name == ref[i,]$gene), ]

                temp = temp[which(str_detect(temp$ref_aa,
                                             aa_variants[[i]][1] )), ]

                temp = temp[which(str_detect(temp$alt_aa,
                                             aa_variants[[i]][2])), ]

                who = which(temp$codon_num == codon_num[[i]] |
                                temp$codon_num == codon_num[[i]] + 1 |
                                temp$codon_num == codon_num[[i]] - 1)

                temp = temp[who, ]

            } else {

                temp = vcf_table[which(vcf_table$Gene_Name == ref[i,]$gene), ]

                temp = temp[which(str_detect(temp$AA_alt, "del")), ]

                who = which(temp$codon_start == codon_num[[i]][1] |
                                temp$codon_start == codon_num[[i]][2])

                temp = temp[who, ]

                if(nrow(temp) != 0) {

                    temp$AA_alt = ref[i, ]$`amino acid`

                }

            }


            voc_data[[ref[i,]$barcode]] = temp



        }

        not_overlapping_variants = names(which(lapply(voc_data, nrow) == 0))


        voc_data = rbindlist(voc_data)

        voc_data = unique(voc_data)



        # collapse table -------------------------------------------------------

        voc_data = voc_data[, c("POS",
                                "DP",
                                "AD_alt",
                                "Gene_Name",
                                "AA_alt",
                                "sample"), with = FALSE]

        voc_data = voc_data[, .(DP = unique(DP),
                                AD_alt = sum(AD_alt)),
                            by = .(POS, Gene_Name, AA_alt, sample)]

        voc_data = voc_data[, .(DP = sum(DP),
                                AD_alt = sum(AD_alt)),
                            by = .(Gene_Name, AA_alt, sample)]


        voc_data$AF = voc_data$AD_alt / voc_data$DP

        # add non overlapping rules --------------------------------------------


        not_overlapping_variants = str_split(not_overlapping_variants,
                                             "\\:", simplify = TRUE)

        voc_data = rbind(voc_data,
                         data.table(Gene_Name = rep(not_overlapping_variants[,1],
                                                    length(unique(voc_data$sample))),

                                    AA_alt = rep(not_overlapping_variants[,2],
                                                 length(unique(voc_data$sample))),

                                    sample = unique(voc_data$sample),

                                    DP = 0,

                                    AD_alt = 0,

                                    AF = 0))

        VoC_hits_list[[ref_index]] = voc_data

    }

    for(i in names(VoC_hits_list)) {

        strain = str_split(i, "\\/", simplify = TRUE)
        strain = strain[, ncol(strain)]
        strain = str_remove(strain, "\\.txt")

        VoC_hits_list[[i]]$lineage = strain

    }

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
