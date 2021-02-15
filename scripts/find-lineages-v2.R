library(vcfR)
library(data.table)
library(stringr)

rm(list = ls())



bases = c("A", "G", "C", "T", "-")


#---- create reference based rules 

vcf.file = read.vcfR("L1_S22_L001_freebayes.vcf", verbose = FALSE)
fix = data.table::as.data.table(vcf.file@fix)

rm(vcf.file)

fix.ref = fix[which(is.na(fix$ALT)), ]
fix = fix[which(!is.na(fix$ALT)), ]

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

# fix.ref = fix[which(str_detect(fix$ALT, "<*>")), ]

# comp = c(comp, 
#          comp = paste(fix.ref$POS, "=='", fix.ref$REF, "'", sep = ""))


#----- create alt rules

alt = str_split(fix$ALT, ",", simplify = TRUE)

alt = as.data.table(alt)

alt$POS = fix$POS

alt$REF = fix$REF

# alt[which(alt$V2 == "<*>"), ]$V2 = alt[which(alt$V2 == "<*>"), ]$REF

fix = alt

rm(alt)
# fix$ALT = str_split(fix$ALT, ",", simplify = TRUE)[,1]

for(i in 1:nrow(fix)){
  
  ref.len = str_length(fix[i,]$REF)
  alt.len = str_length(fix[i,]$V1)
  # exclude = c()
  
  if((ref.len == 1) & (alt.len == 1)) {
    
    comp = c(comp,
             paste(fix[i,]$POS, "=='", fix[i,]$V1, "'", sep = ""))
    
    comp = c(comp,
             paste(fix[i,]$POS, "!='", bases[!(bases %in% fix[i,]$V1)], "'", sep = ""))
    
    # exclude = c(exclude, fix[i,]$V1)
    
    if(fix[i,]$V2 != "") {
      
      comp = c(comp,
               paste(fix[i,]$POS, "=='", fix[i,]$V2, "'", sep = ""))
      
      comp = c(comp,
               paste(fix[i,]$POS, "!='", bases[!(bases %in% fix[i,]$V2)], "'", sep = ""))
      
      # exclude = c(exclude, fix[i,]$V2)
      
    }
    
    # comp = c(comp,
    #          paste(fix[i,]$POS, "!='", bases[!(bases %in% exclude)], "'", sep = ""))
    
  }
  
  if((ref.len != 1) & (alt.len == 1)) {
    
    ref.split = str_split(fix[i,]$REF, "", simplify = TRUE)
    
    comp = c(comp, 
             paste(fix[i,]$POS, "=='", 
                   fix[i,]$V1, "'", 
                   sep = ""))
    
    comp = c(comp,
             paste(fix[i,]$POS, "!='", bases[!(bases %in% fix[i,]$V1)], "'", sep = ""))
    
    # exclude = c(exclude, fix[i,]$V1)
    
    if(fix[i,]$V2 != "") {
      comp = c(comp, 
               paste(as.numeric(fix[i,]$POS) + 0:(ref.len - 1), 
                     "=='", 
                     ref.split, "'", 
                     sep = ""))
      
      for(j in 1:ref.len) {
        
        comp = c(comp,
                 paste(as.numeric(fix[i,]$POS) + j - 1, "!='", 
                       bases[!(bases %in% ref.split[j])], "'", 
                       sep = ""))
        
      }
      
      # exclude = c(exclude, ref.split)
      
    }
    
    # comp = c(comp,
    #          paste(fix[i,]$POS, "!='", bases[!(bases %in% exclude)], "'", sep = ""))
    
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

rm(fix, fix.ref, ref.split, temp, alt.len, bases, i, ref.len, one.run)

#---- read decision rules

decisions = data.table::fread("decision_tree_rules_2021_01_13.txt", header = FALSE)

rules = str_split(decisions$V2, ",", simplify = TRUE)
rules = as.data.table(rules)

one.run = function(x) {
  
  return(length(which(x != "")))
  
}

rules$n.rules = apply(rules, 1, one.run)

out = list()
k = 1

for(i in 1:(ncol(rules) - 1)) {

  who = rules[[i]] == ""
  
  if(length(which(who)) > 0){
    
    temp = decisions[which(who), ]
    temp$overlap = i - 1
    temp$total = rules[which(who), ]$n.rules
    
    temp$ratio = temp$overlap / temp$total
    
    out[[k]] = temp
    
    k = k + 1
    
    decisions = decisions[which(!who), ]
    rules = rules[which(!who), ]
    
  }
  
  if(nrow(decisions) > 0) {
    who = rules[[i]] %in% comp$whole.expr
    
    temp = decisions[which(!who), ]
    
    if(nrow(temp) > 0) {
      temp$overlap = i - 1
      temp$total = rules[which(!who), ]$n.rules
      
      temp$ratio = temp$overlap / temp$total
      
      out[[k]] = temp
      
      k = k + 1
      
    }

    decisions = decisions[which(who), ]
    rules = rules[which(who), ]
    
    cat(c("Rule:", i, "\n"))
  } else {
    
    break
  
  }
  
}


out = rbindlist(out)

out = out[order(out$ratio, decreasing = TRUE), ]

colnames(out) = c("Lineage", "Rules", "Overlap", "Total", "Ratio")

write.table(out, "overlapping-lineages-v4.tsv", row.names = FALSE, quote = FALSE, sep = "\t")
