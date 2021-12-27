#' merge_vcf
#'
#' @description
#' Merge Variant Calling Format (VCF) files into a single tab-delimited table
#'
#' @param vcf_fls
#' A list of paths to VCF files
#'
#' @param vcf_folder
#' A path to a folder containing all VCF file that
#' will be integrated into a single table
#'
#' @param print.out
#' Logical value indicating if the produced table should be printed
#'
#' @param file.out
#' Given name for the output file
#'
#' @param gff3_path
#' Path to GFF3 file
#'
#'
#' @return
#' A data table contaiing all variants from each sample of the input VCF files
#'
#' @export
#'
#' @examples
#'
#' merge_vcf(vcf_folder = system.file("extdata",
#'                                    "vcf-files",
#'                                    package = "lineagespot"),
#'
#'           gff3_path = system.file("extdata",
#'                                   "NC_045512.2_annot.gff3",
#'                                   package = "lineagespot"))


merge_vcf <- function(vcf_fls = NULL,
                      vcf_folder = NULL,
                      gff3_path = NULL,
                      file.out = paste0("Variant_table_", Sys.Date(), ".txt"),
                      print.out = FALSE) {


    if( is.null(vcf_fls) & is.null(vcf_folder) ) {

        stop("Please provide some VCF files")

    }


    if( is.null(vcf_fls) ) {

        vcf_fls = list.files(vcf_folder,
        pattern = "vcf",
        full.names = TRUE)


    }

    if( is.null(gff3_path) ){

        stop("Please provide a valid GFF file containing gene coordinates")

    }


    # Read input VCF files -----------------------------------------------------

    vcf_list = list()

    for(i in vcf_fls) {

        # vcf_list[[i]] = vcfR::read.vcfR(i, verbose = FALSE)
        
        vcf_list[[i]] = readVcf(i)

    }

    rm(i)


    # Convert VCF to table -----------------------------------------------------

    vcf_list = lapply(vcf_list, vcf_to_table)

    # break multiple rules -----------------------------------------------------

    vcf_list = lapply(vcf_list, break_multiple_variants)

    # add variant parameters ---------------------------------------------------

    vcf_list = lapply(vcf_list, compute_AF)

    # add variant parameters ---------------------------------------------------

    vcf_list = lapply(vcf_list, change_AA_abbreviations)

    # add variant parameters ---------------------------------------------------

    genes = read_gene_coordinates(gff_path = gff3_path)

    vcf_list =lapply(vcf_list, correct_Orf1ab_gene, genes)

    # add sample name ----------------------------------------------------------

    for(i in names(vcf_list)) {

        sample_name = unlist(str_split(i, "\\/"))
        sample_name = sample_name[length(sample_name)]

        sample_name = str_split(sample_name, "\\.vcf", simplify = TRUE)[,1]

        vcf_list[[i]]$sample = sample_name

    }

    # Combine all vcfs into one table ------------------------------------------

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

    # out = cbind(x@fix[, seq_len(7)],
    # vcfR::extract_gt_tidy(x, verbose = FALSE),
    # vcfR::extract_info_tidy(x))
    
    out = cbind(as.data.frame(rowRanges(x)),
                as.data.frame(fixed(x)),
                as.data.frame(info(x)))
    
    out$ALT = unlist(lapply(out$ALT, function(x){ return(paste(as.character(x), collapse = ",")) }))
    
    out$ANN = unlist(lapply(out$ANN, function(x){ return(x[1]) }))
    
    out$ID = row.names(out)
    
    out$AD = unlist(lapply(assays(x)[["AD"]], function(x) {return(paste(x, collapse = ","))}))

    out = setDT(out)

    out = out[which(out$ALT != "<*>"), ]

    ANN_matrix = str_split(out$ANN, "\\|", simplify = TRUE)[,c(4, 10, 11)]

    # ANN_matrix = paste0(ANN_matrix[,1], "|", ANN_matrix[,2])

    values_wewant = c("seqnames", # "CHROM",
                      "start", # "POS",
                      "ID",
                      "REF",
                      "ALT",
                      "DP",  # "gt_DP",
                      "AD" ) # "gt_AD")

    result = intersect(values_wewant, colnames(out))

    out = out[, result, with = FALSE]

    if(identical((intersect(result, "gt_DP")), "gt_DP")){

        names(out)[names(out) == "gt_DP"] = "DP"

    }

    if(identical((intersect(result, "gt_AD")), "gt_AD")){

        names(out)[names(out) == "gt_AD"] = "AD"

    }

    colnames(out) = c("CHROM", "POS", "ID", "REF", "ALT", "DP", "AD")

    out$CHROM = str_squish(out$CHROM)
    out$Gene_Name = ANN_matrix[,1]
    out$Nt_alt = ANN_matrix[,2]
    out$AA_alt = ANN_matrix[,3]

    return(out)

}

break_multiple_variants <- function(x) {

    ALT_matrix = str_split(x$ALT, ",", simplify = TRUE)
    AD_matrix = str_split(x$AD, ",", simplify = TRUE)

    out = list()

    for(i in seq_len(ncol(ALT_matrix))) {

        who = which(ALT_matrix[,i] != "")


        out[[i]] = data.table(CHROM = x[who, ]$CHROM,
                              POS = x[who, ]$POS,
                              ID = x[who, ]$ID,
                              REF = x[who, ]$REF,
                              ALT = ALT_matrix[who, i],
                              DP = x[who, ]$DP,
                              AD_ref = AD_matrix[who, 1],
                              AD_alt = AD_matrix[who, i+1],
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

    for(i in seq_len(nrow(AA_abbreviations))) {

        x$AA_alt = str_replace_all(x$AA_alt,
        AA_abbreviations[i,]$Three_Letter,
        AA_abbreviations[i,]$One_Letter)

    }

    return(x)

}

correct_Orf1ab_gene <- function(x, genes) {

    x$codon_num = 0

    x = x[order(x$POS), ]

    for(i in seq_len(nrow(genes))) {

        who = which(x$POS >= genes[i,]$start_pos & x$POS <= genes[i,]$end_pos)

        x[who, ]$codon_num = x[who, ]$POS - genes[i,]$start_pos + 1
        x[who, ]$Gene_Name = genes[i,]$gene_name
        x[who, ]$codon_num = ( x[who, ]$codon_num %/% 3 ) +
                             ( x[who, ]$codon_num %% 3 > 0 )

    }

    x = x[which(x$codon_num != 0), ]

    return(x)

}

read_gene_coordinates <- function(gff_path) {

    gene_annot = fread(gff_path, header = FALSE, sep = "\t")

    gene_annot = gene_annot[, c(4, 5, 9), with = FALSE]

    colnames(gene_annot) = c("start_pos", "end_pos", "gene_name")

    gene_annot[which(gene_annot$gene_name %in% paste0("NSP", seq_len(10))), ]$gene_name = "ORF1a"

    gene_annot[which(str_detect(gene_annot$gene_name, "NSP")), ]$gene_name = "ORF1b"

    gene_annot = gene_annot[, .(start_pos = min(start_pos),
    end_pos = max(end_pos)), by = gene_name]

    gene_annot = gene_annot[order(gene_annot$start_pos), ]

    return(gene_annot)


}
