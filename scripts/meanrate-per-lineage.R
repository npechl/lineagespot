library(data.table)
library(stringr)

library(dplyr)

rules.comp = fread("overlapping-lineages.tsv")

lineage.average = rules.comp %>% group_by(Lineage) %>% summarise(n = mean(Overlap.rate))

lineage.average = lineage.average[order(lineage.average$n, decreasing = TRUE), ]

write.table(lineage.average, "mean-per-lineage.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

