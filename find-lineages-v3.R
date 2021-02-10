library(vcfR)
library(data.table)
library(stringr)
library(seqinr)

# library(stringdist)

rm(list = ls())


bases = c("A", "G", "C", "T", "-")


#---- create reference based rules 

ref = seqinr::read.fasta("NC_045512.fasta")
ref = as.character(ref[[1]])

ref = data.table(pos = 1:length(ref),
                 logic.expr = "==",
                 base = str_to_upper(ref))

ref$whole.expr = paste(ref$pos, ref$logic.expr, "'", ref$base, "'", sep = "")

#---- find variances

vcf.file = read.vcfR("vcf-files/L1_S22_L001_sars_variants_filtered.vcf", verbose = FALSE)
fix = data.table::as.data.table(vcf.file@fix)

rm(vcf.file)

fix = fix[which(fix$ALT != "<*>"), ]
fix$ALT = str_remove_all(fix$ALT, "\\,<\\*>")

ALT = str_split(fix$ALT, ",", simplify = TRUE)

out = list()

for(i in 1:ncol(ALT)) {
  
  who = which(ALT[,i] != "")
  
  out[[i]] = data.table(pos = fix[who, ]$POS,
                        ref = fix[who, ]$REF,
                        alt = ALT[who, i])
  
}

fix = rbindlist(out)

rm(out, ALT, i, who)

comp = list()

exclude = c()

for(i in 1:nrow(fix)){
  
  ref.len = str_length(fix[i,]$ref)
  alt.len = str_length(fix[i,]$alt)
  # exclude = c()
  
  if((ref.len == 1) & (alt.len == 1)) {
    
    temp = data.table(pos = fix[i,]$pos,
                      logic.expr = "==",
                      base = fix[i,]$alt)
    
    temp$whole.expr = paste(temp$pos, temp$logic.expr, "'", temp$base, "'", sep = "")
    
    comp[[i]] = temp
    
  } else if((ref.len != 1) & (alt.len == 1)) {
    
    ref.split = str_split(fix[i,]$ref, "", simplify = TRUE)
    
    temp = data.table(pos = c(as.numeric(fix[i,]$pos), as.numeric(fix[i,]$pos) + 1:(ref.len - 1)),
                      logic.expr = "==",
                      base = c(fix[i,]$alt, rep("-", ref.len - 1)))
    
    temp$whole.expr = paste(temp$pos, temp$logic.expr, "'", temp$base, "'", sep = "")
    
    comp[[i]] = temp
    
  } else if((ref.len != 1) & (alt.len != 1) & (ref.len == alt.len)) {
    
    ref.split = str_split(fix[i,]$ref, "", simplify = TRUE)
    alt.split = str_split(fix[i,]$alt, "", simplify = TRUE)
    
    ref.split = as.vector(ref.split)
    alt.split = as.vector(alt.split)
    
    diff = stringdist::stringdist(ref.split, alt.split)
    
    who = which(diff != 0)
    
    temp = data.table(pos = c(as.numeric(fix[i,]$pos) + who - 1),
                      logic.expr = "==",
                      base = alt.split[who])
    
    temp$whole.expr = paste(temp$pos, temp$logic.expr, "'", temp$base, "'", sep = "")
    
    comp[[i]] = temp
    
  } else if((ref.len != 1) & (alt.len != 1) & (ref.len > alt.len)) {
    
    ref.split = str_split(fix[i,]$ref, "", simplify = TRUE)
    alt.split = str_split(fix[i,]$alt, "", simplify = TRUE)
    
    ref.split = as.vector(ref.split)
    alt.split = as.vector(alt.split)
    
    alt.j = 1
    
    gaps = c()
    
    for(j in 1:length(ref.split)) {
      
      if(ref.split[j] != alt.split[alt.j]) {
        
        gaps = c(gaps, 
                 as.numeric(fix[i,]$pos) + j - 1)
        
      } else {
        
        alt.j = ifelse((alt.j + 1) > alt.len, alt.len, alt.j + 1)
        
      }
      
    }
    
    if(is.null(gaps)) {
      
      gaps = as.numeric(fix[i,]$pos) + 1
      
    }
    
    temp = data.table(pos = gaps,
                      logic.expr = "==",
                      base = "-")
    
    temp$whole.expr = paste(temp$pos, temp$logic.expr, "'", temp$base, "'", sep = "")
    
    comp[[i]] = temp
    
    
  } else if((ref.len != 1) & (alt.len != 1) & (ref.len < alt.len)) {
    
    ref.split = str_split(fix[i,]$ref, "", simplify = TRUE)
    alt.split = str_split(fix[i,]$alt, "", simplify = TRUE)
    
    ref.split = as.vector(ref.split)
    alt.split = as.vector(alt.split)
    
    ref.j = ref.len
    
    insertions = c()
    
    for(j in length(alt.split):1) {
      
      if(alt.split[j] != ref.split[ref.j]) {
        
        insertions = c(insertions, j)
        
      } else {
        
        ref.j = ifelse((ref.j - 1) < 0, 1, ref.j - 1)
        
      }
      
    }
    
    if(length(insertions) == 1) {
      
      temp = data.table(pos = c(as.numeric(fix[i,]$pos) + insertions - 1),
                        logic.expr = "==",
                        base = alt.split[insertions])
      
      temp$whole.expr = paste(temp$pos, temp$logic.expr, "'", temp$base, "'", sep = "")
      
      comp[[i]] = temp
      
    } else {
      
      exclude = c(exclude, i)
    }
  
  } else {
    
    exclude = c(exclude, i)
    
  }
  
}

comp = rbindlist(comp)

fix = comp

rm(list=setdiff(ls(), c("fix", "ref", "bases")))

fix = rbind(ref[which(!(ref$pos %in% fix$pos)), ],
            fix)

fix$pos = as.numeric(fix$pos)

fix = fix[order(fix$pos), ]


one.run = function(x, bases) {
  
  # all = all[which(all$pos == x), ] 
  
  out = data.table(pos = as.numeric(x[1]),
                   logic.expr = "!=",
                   base = bases[!(bases %in% x[3])])
  
  out$whole.expr = paste(out$pos, out$logic.expr, "'", out$base, "'", sep = "")
  
  return(out)
}

# pos = unique(fix$pos)

comp = apply(fix, 1, one.run, bases = bases)

comp = rbindlist(comp)

comp = rbind(fix, comp)

#---- read decision rules

rm(bases, one.run, fix, ref)

decisions = data.table::fread("decision_tree_rules_2021_01_13.txt", header = FALSE)

rules = str_split(decisions$V2, ",", simplify = TRUE)
rules = as.data.table(rules)

one.run = function(x, expr) {
  
  x = x[which(x != "")]
  
  out = data.table(total.len = length(x),
                   total.overlap = length(which(expr %in% x)))
  
  return(out)
  
}

rules.stats = apply(rules, 1, one.run, comp$whole.expr)
rules.stats = rbindlist(rules.stats)

rules$n.rules = rules.stats$total.len

rules$n.totalOverlap = rules.stats$total.overlap

out = list()
k = 1

for(i in 1:(ncol(rules) - 2)) {
  
  who = rules[[i]] == ""
  
  if(length(which(who)) > 0){
    
    temp = decisions[which(who), ]
    temp$total = rules[which(who), ]$n.rules
    temp$tree.overlap = i - 1
    temp$total.overlap = rules[which(who), ]$n.totalOverlap
    
    # temp$ratio = temp$overlap / temp$total
    
    out[[k]] = temp
    
    k = k + 1
    
    decisions = decisions[which(!who), ]
    rules = rules[which(!who), ]
    
  }
  
  if(nrow(decisions) > 0) {
    who = rules[[i]] %in% comp$whole.expr
    
    temp = decisions[which(!who), ]
    
    if(nrow(temp) > 0) {
      temp$total = rules[which(!who), ]$n.rules
      temp$tree.overlap = i - 1
      temp$total.overlap = rules[which(!who), ]$n.totalOverlap
      
      # temp$ratio = temp$overlap / temp$total
      
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

out$tree.ratio = out$tree.overlap / out$total
out$total.ratio = out$total.overlap / out$total

out = out[order(out$tree.ratio, decreasing = TRUE), ]

colnames(out) = c("Lineage", "Rules", "Total", "Tree.Overlap", "Total.Overlap", "Tree.Ratio", "Total.Ratio")

write.table(out, "OverLine-L1_S22_L001_sars_variants_filtered.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

rm(list=setdiff(ls(), c("out", "comp", "decisions")))
