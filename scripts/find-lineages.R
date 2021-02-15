library(vcfR)
library(data.table)
library(stringr)

rm(list = ls())



bases = c("A", "G", "C", "T", "-")



vcf.file = read.vcfR("CoV19_L1_S1_sorted_uniq_mpileup.vcf", verbose = FALSE)
fix = data.table::as.data.table(vcf.file@fix)

rm(vcf.file)

fix.ref = fix[which(fix$ALT == "<*>"), ]
fix = fix[which(fix$ALT != "<*>"), ]

comp = paste(fix.ref$POS, "=='", fix.ref$REF, "'", sep = "")


one.run = function(x, bases) {
  
  x.split = str_split(x, "='", simplify = TRUE)
  pos = str_sub(x.split[,1], 1, -2)
  curr.base = str_sub(x.split[,2], 1, -2)
  
  out = c(x, paste(pos, "!='", bases[!(bases %in% curr.base)], "'", sep = ""))
  
  return(out)
}

comp = lapply(as.list(comp), one.run, bases = bases)
comp = unlist(comp)

fix$ALT = str_split(fix$ALT, ",", simplify = TRUE)[,1]

for(i in 1:nrow(fix)){
  
  ref.len = str_length(fix[i,]$REF)
  alt.len = str_length(fix[i,]$ALT)
  
  if((ref.len == 1) & (alt.len == 1)) {
    
    comp = c(comp,
             paste(fix[i,]$POS, "=='", fix[i,]$ALT, "'", sep = ""))
    
    comp = c(comp,
             paste(fix[i,]$POS, "!='", bases[!(bases %in% fix[i,]$ALT)], "'", sep = ""))
    
  }
  
  if((ref.len != 1) & (alt.len == 1)) {
    
    ref.split = str_split(fix[i,]$REF, "", simplify = TRUE)
    
    comp = c(comp, 
             paste(fix[i,]$POS, "=='", 
                   ref.split[1], "'", sep = ""))
    
    comp = c(comp,
             paste(fix[i,]$POS, "!='", bases[!(bases %in% fix[i,]$ALT)], "'", sep = ""))
    
    comp = c(comp, 
             paste(as.numeric(fix[i,]$POS) + 1:(ref.len - 1), "=='-'", sep = ""))
    
    
  }
  
}






comp = str_split(comp, ",")
comp = unlist(comp)
comp = unique(comp)

temp = str_locate(comp, "='")

comp = data.table(pos = as.numeric(str_sub(comp, 1, temp[,1] - 2)),
                  logix.expr = str_sub(comp, temp[,1] - 1, temp[,2] - 1),
                  base = str_sub(comp, temp[,2], -1))

comp$whole.expr = paste(comp$pos, 
                        comp$logix.expr, 
                        comp$base, 
                        sep = "")



decisions = data.table::fread("decision_tree_rules_2021_01_13_sorted.txt")

out = list()

for(i in unique(decisions$lineage)) {
  
  temp = decisions[which(decisions$lineage == i), ]
  
  rules.temp = str_split(temp$rules, ",")
  rules.temp = unlist(rules.temp)
  rules.temp = unique(rules.temp)
  
  out[[i]] = data.table(lineage = i,
                        rules = paste(rules.temp, collapse = ","))
  
}

out = rbindlist(out)
decisions = out
rm(out)

rules = decisions$rules
rules = str_split(rules, ",")
rules = unlist(rules)
rules = unique(rules)

temp = str_locate(rules, "='")

rules = data.table(pos = as.numeric(str_sub(rules, 1, temp[,1] - 2)),
                   logix.expr = str_sub(rules, temp[,1] - 1, temp[,2] - 1),
                   base = str_sub(rules, temp[,2], -1))

rm(temp)

rules = rules[order(rules$pos), ]

rules$whole.expr = paste(rules$pos, 
                         rules$logix.expr, 
                         rules$base, 
                         sep = "")

comp = comp$whole.expr




# rules.over = list()

one.run = function(x, comp) {
  
  x.split = str_split(x[2], ",", simplify = TRUE)
  
  comp.over = comp[which(comp %in% x.split)]
  
  nover = length(comp.over)
  
  if(nover == 0){
    comp.over = ""
  } else {
    comp.over = paste(comp.over, collapse = ";")
  }
  
  x.split = t(x.split)
  
  # return(x.split)
  
  return(data.table(Lineage = x[1],
                    Overlap = nover,
                    Lineage.rules.number = length(x.split),
                    Overlap.rate = (nover / length(x.split)),
                    Lineage.rules = x[2],
                    Overlapping.rules = comp.over))

}


rules.over = apply(decisions, 1, one.run, comp)

rules.over = data.table::rbindlist(rules.over)

rules.over = unique(rules.over)

rules.over = rules.over[which(rules.over$Overlap != 0), ]
rules.over = rules.over[order(rules.over$Lineage, rules.over$Overlap.rate, decreasing = TRUE), ]

rules.over$Lineage.rules = str_replace_all(rules.over$Lineage.rules, ",", ";")


write.table(rules.over, "overlapping-lineages-v2.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

rm(list = ls())

rules.over = fread("overlapping-lineages.csv")

rules = data.table(Lineage = rules.over$Lineage,
                   Rule.index = paste("Rule", 1:nrow(rules.over), sep = ""),
                   Lineage.rules = rules.over$Lineage.rules)


rules.comp = combn(rules$Rule.index, 2)
rules.comp = t(rules.comp)
rules.comp = as.data.table(rules.comp)

who = match(rules.comp$V1, rules$Rule.index)
rules.comp$V3 = rules[who,]$Lineage.rules

who = match(rules.comp$V2, rules$Rule.index)
rules.comp$V4 = rules[who,]$Lineage.rules

rm(rules.over, who)

one.run = function(x) {
  
  x.split.1 = str_split(x[3], ";", simplify = TRUE)
  x.split.2 = str_split(x[4], ";", simplify = TRUE)
  
  return(length(which(x.split.1 %in% x.split.2)))
  
}

nover.rules = apply(rules.comp, 1, one.run)

one.run = function(x) {
  
  x.split.1 = str_split(x[3], ";", simplify = TRUE)
  
  return(length(x.split.1))
  
}

rules$number.of.rules = apply(rules, 1, one.run)

write.table(rules, "lineage-index.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

who = match(rules.comp$RuleV1, rules$Rule.index)
rules.comp$RuleV1.total = rules[who, ]$number.of.rules

who = match(rules.comp$RuleV2, rules$Rule.index)
rules.comp$RuleV2.total = rules[who, ]$number.of.rules

write.table(rules.comp, "lineage-rules-comparison-v2.tsv", row.names = FALSE, quote = FALSE, sep = "\t")




