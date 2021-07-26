

#' Lineage report for variants overlapping
#'
#' @param hits_table
#'
#' A tab-delimited table containing the overlaps identified
#'
#' @param print.out
#'
#' logical value indicating if the produced table should be printed
#'
#' @param AF_threshold
#'
#' A parameter indicating the AF threshold for identifying variants per sample
#'
#' @import data.table
#' @import stringr
#'
#' @return
#' @export uniq_variants
#'
#' @examples


uniq_variants <- function(hits_table = NULL,
                          AF_threshold = 0.8,
                          print.out = FALSE) {

  if( base::is.null(hits_table) ) {

    base::stop('Please provide a tab-delimited table containing SARS-CoV-2 overlaps')

  }

  # lineage overview --------------------------------------------------------------

  lineage_stats = unique(hits_table[, c("Gene_Name", "AA_alt", "lineage"), with = FALSE])
  lineage_stats = lineage_stats[, .N, by = lineage]

  count_lineages = hits_table[which(hits_table$AF >= AF_threshold), ]

  count_lineages = unique(count_lineages[, c("Gene_Name", "AA_alt", "sample", "lineage"), with = FALSE])

  count_lineages = count_lineages[, .N, by = .(lineage, sample)]

  # mean AF of all variants -------------------------------------------------------

  meanAF = hits_table[, .(meanAF = mean(AF)), by = .(lineage, sample)]


  # Find unique variants per lineage ----------------------------------------------

  count_uniq = hits_table[, c("Gene_Name", "AA_alt", "lineage"), with = FALSE]
  count_uniq = unique(count_uniq)
  count_uniq = count_uniq[order(count_uniq$Gene_Name, count_uniq$AA_alt), ]
  count_uniq = count_uniq[, .N, by = .(Gene_Name, AA_alt)]


  count_uniq = count_uniq[which(count_uniq$N == 1), ]

  hits_table_uniq = hits_table[which(hits_table$Gene_Name %in% count_uniq$Gene_Name & hits_table$AA_alt %in% count_uniq$AA_alt), ]

  # mean AF of unique variants --------------------------------------------------------------

  meanAF_uniq = hits_table_uniq[, .(meanAF_uniq = mean(AF)), by = .(lineage, sample)]

  # filter out zero-AF variants -----------------------------------------------------

  hits_table_uniq = hits_table_uniq[which(hits_table_uniq$AF != 0), ]

  # min AF of non zeror unique variants

  minAF_uniq = hits_table_uniq[, .(minAF_uniq_nonzero = min(AF)), by = .(lineage, sample)]

  # create output table -------------------------------------------------------------

  overall = rbindlist(list(meanAF, meanAF_uniq, minAF_uniq, count_lineages),
                      use.names = TRUE,
                      fill = TRUE)

  overall = overall[order(overall$lineage, overall$sample), ]

  overall = overall[, lapply(.SD, sum, na.rm = TRUE), by = .(lineage, sample)]

  overall[which(overall$minAF_uniq_nonzero == 0), ]$minAF_uniq_nonzero = NA

  who = base::match(overall$lineage, lineage_stats$lineage)

  overall$`lineage N. rules` = lineage_stats[who, ]$N

  overall$`lineage prop.` = overall$N / overall$`lineage N. rules`

  if( print.out ) {

    data.table::fwrite(overall,
                       file = "lineage-report.txt",
                       row.names = FALSE,
                       quote = FALSE,
                       sep = "\t")

  }

  return(overall)



}


