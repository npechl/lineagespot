library(xlsx)
library(stringr)
library(data.table)

rm(list = ls())

source("scripts/help_functions.R")

run_2021_02_01 = read.xlsx("human/Analysis-COVID19-run-2021-02-01.xlsx", sheetIndex = 1)
lin.names = run_2021_02_01[3, ]
run_2021_02_01 = run_2021_02_01[4:nrow(run_2021_02_01),]
colnames(run_2021_02_01) = lin.names
run_2021_02_01 = as.data.table(run_2021_02_01)

run_2021_02_05 = read.xlsx("human/Analysis-COVID19-run-2021-02-05.xlsx", sheetIndex = 1)
lin.names = run_2021_02_05[3, ]
run_2021_02_05 = run_2021_02_05[4:nrow(run_2021_02_05),]
colnames(run_2021_02_05) = lin.names
run_2021_02_05 = as.data.table(run_2021_02_05)

rm(lin.names)


sewage.lin = fread("vcf-files/Sewage-L2_S10_L001_freebayes-lineagespot.tsv")
# 
# sewage.lin.f = sewage.lin[which(sewage.lin$Total.Ratio >= 0.75), ]
# sewage.lin.f = sewage.lin.f[which(sewage.lin.f$Total.Avg.AD != 0), ]

sewage.lin.f = sewage.lin

sewage.lin.f.collapsed = collapse_lineages_table(sewage.lin.f)

cont_run_1 = table(run_2021_02_01$Lineage)
cont_run_1 = sort(cont_run_1, decreasing = TRUE)
cont_run_1 = as.data.frame(cont_run_1)
cont_run_1$run = "run-2021-02-01"
colnames(cont_run_1) = c("Lineage", "Freq", "Human run")

who = match(cont_run_1$Lineage, sewage.lin.f.collapsed$Lineage)

cont_run_1 = cbind(cont_run_1, sewage.lin.f.collapsed[who, 2:ncol(sewage.lin.f.collapsed)])


cont_run_2 = table(run_2021_02_05$Lineage)
cont_run_2 = sort(cont_run_2, decreasing = TRUE)
cont_run_2 = as.data.frame(cont_run_2)
cont_run_2$run = "run-2021-02-05"
colnames(cont_run_2) = c("Lineage", "Freq", "Human run")

who = match(cont_run_2$Lineage, sewage.lin.f.collapsed$Lineage)

cont_run_2 = cbind(cont_run_2, sewage.lin.f.collapsed[who, 2:ncol(sewage.lin.f.collapsed)])

cont_run = rbind(cont_run_1, cont_run_2)

xlsx::write.xlsx(cont_run, file = "comparisons.xlsx", sheetName = "human-compare", row.names = FALSE, showNA = FALSE)


old_mink = fread("vcf-files/L1_S22_L001_freebayes-lineagespot-new.tsv")

old_mink = old_mink[which(old_mink$Total.Ratio >= 0.75), ]
old_mink = old_mink[which(old_mink$Total.Avg.AD != 0), ]

old_mink = collapse_lineages_table(old_mink)

who = match(old_mink$Lineage, sewage.lin.f.collapsed$Lineage)

old_mink = cbind(old_mink[,1], sewage.lin.f.collapsed[who, 2:ncol(sewage.lin.f.collapsed)])

xlsx::write.xlsx(old_mink, file = "comparisons.xlsx", sheetName = "mink-compare", row.names = FALSE, showNA = FALSE, append = TRUE)

