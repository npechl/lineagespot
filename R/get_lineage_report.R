#' get_lineage_report
#'
#' Retrieving lineage varinats related information from outbreak.info API
#'
#' @param lineages
#' a character vector containing the names of the lineages of interest
#'
#' @param base.url
#' The base API URL used to search for lineage reports
#'
#'
#' @return A list of data table elements of lineage reports
#' @export get_lineage_report
#'
#' @examples
#'
#' get_lineage_report(lineages = c("B.1.1.7", "B.1.617.2"))
#'
#'



get_lineage_report <- function(lineages, base.url = "https://api.outbreak.info/genomics/lineage-mutations?pangolin_lineage=") {

    lineages = str_to_upper(lineages)

    out = list()

    for(l in lineages) {

        strain = GET(paste0(base.url, l))

        strain = content(strain)

        if(strain$success) {

            strain = strain$results[[l]]
            strain = rbindlist(strain)

            strain$`amino acid` = str_split(strain$mutation, "\\:",
                                            simplify = TRUE)[,2]

            strain[which(strain$type != "deletion"), ]$`amino acid` = str_to_upper(strain[which(strain$type != "deletion"), ]$`amino acid`)

            strain$barcode = paste0(strain$gene, ":", strain$`amino acid`)

            strain[which(strain$type != "deletion"), ]$type = "snp"

            out[[l]] = strain

        } else {

            message(c("No data has been retrieved for ", l))
            return(NULL)

        }

    }


    return(out)




}
