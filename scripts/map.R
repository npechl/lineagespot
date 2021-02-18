# clear
cat("\014")
rm(list = ls())
dev.off(dev.list()["RStudioGD"])

# source
source("help_functions.R")

# Libraries
library(openxlsx)

# INFO
#INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus"
#FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele"

# ---------------------- PARAMETERS --------------------------------------------
# Specify thresholds for tree ration and total ratio
total_ratio_threshold <- 0.75
tree_ratio_threshold <- 0 # We keep the rows that have Total.Tree.Ratio >= tree_ratio_threshold (meaningless here)

#-------------------- INPUT FILES ----------------------------------------------
samples_detected_filepath <- 'lineages_detected/INFO FOR FOTIS.xlsx' # (fotis file)

analysis_output_folder <- 'new-vcf-files'
analysis_output_filename <- 'L1_S22_L001_freebayes-lineagespot-new.tsv'
analysis_output_filepath <- paste(analysis_output_folder, analysis_output_filename,sep = '/')


#---------------------- Importing detected lineages ----------------------------
# Detected lineages
lineages_detected <- read.xlsx(samples_detected_filepath, sheet = 1, startRow = 1, colNames = TRUE, rowNames = FALSE)
lineages_detected$DATE.OF.RECEIPT <- as.Date(lineages_detected$DATE.OF.RECEIPT, format = "%d.%m.%y")
lineages_detected <- lineages_detected[order(lineages_detected$DATE.OF.RECEIPT), ] # Sort by date

# Separate by months
lineages_october <- subset(lineages_detected, format.Date(DATE.OF.RECEIPT, "%m")=="10")
lineages_november <- subset(lineages_detected, format.Date(DATE.OF.RECEIPT, "%m")=="11")
lineages_december <- subset(lineages_detected, format.Date(DATE.OF.RECEIPT, "%m")=="12")

# count percentages of occurences
october_lineages_perc <- perc_of_lineages(lineages_october$Lineage)
november_lineages_perc <- perc_of_lineages(lineages_november$Lineage)
december_lineages_perc <- perc_of_lineages(lineages_december$Lineage)


# ------------------ Importing results from our analysis -----------------------

# Import lineages from our analysis
filename <- strsplit(analysis_output_filename, '.', fixed = TRUE)[[1]][1]
lineages_pred <- read.csv(analysis_output_filepath, sep = '\t')

# Filter lineages,keep only those that have tree ratio >= 0 (absolutely nothing :P)
lineages_pred <- lineages_pred[which(lineages_pred$Total.Avg.AD != 0), ]
lineages_pred <- lineages_pred[which(lineages_pred$Tree.Ratio >= tree_ratio_threshold), ]

# Collapse
lineages_collapsed <- collapse_lineages_table(lineages_pred)

# Add new columns to initial table
lineages_pred <- add_columns_to_initial_matrix(lineages_pred, lineages_collapsed )

# Save tables
utils::write.table(lineages_pred, file = paste(analysis_output_folder, '/', filename, '_modified.tsv', sep = ''), 
                   row.names = FALSE, 
                   quote = FALSE, 
                   sep = "\t")

utils::write.table(lineages_collapsed, file = paste(analysis_output_folder, '/', filename, '_collapsed.tsv', sep = ''), 
                   row.names = FALSE, 
                   quote = FALSE, 
                   sep = "\t")