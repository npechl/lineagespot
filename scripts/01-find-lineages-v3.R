# Libraries needed ----------------------------------

# library(vcfR)
# library(data.table)
# library(stringr)
# library(seqinr)
# library(stringdist)

base::rm(list = ls())


# Inputs ----------------------------------

bases = base::c("A", "G", "C", "T", "-")

reference.path = "ref/NC_045512.fasta"

vcf.path = "vcf-files/Sewage_CoV19_L3_S1_s500_freebayes.vcf"

decision.rules.path = "ref/decision_tree_rules.txt"

nreads = 490549


utils::download.file(
  url = "https://raw.githubusercontent.com/cov-lineages/pangoLEARN/master/pangoLEARN/data/decision_tree_rules.txt",
  destfile = decision.rules.path,
  method = "curl"
)


# Create reference based rules ----------------------------------

start.time = base::Sys.time()

ref = seqinr::read.fasta(reference.path)
ref = base::as.character(ref[[1]])

ref = data.table::data.table(pos = 1:base::length(ref),
                             logic.expr = "==",
                             base = stringr::str_to_upper(ref))

ref$whole.expr = base::paste(ref$pos, 
                             ref$logic.expr, 
                             "'", ref$base, "'", 
                             sep = "")

# Create variant based rules ----------------------------------

vcf.file = vcfR::read.vcfR(vcf.path, verbose = FALSE)
fix = data.table::as.data.table(vcf.file@fix)

gt.tidy = vcfR::extract_gt_tidy(vcf.file, verbose = FALSE)

gt.ad = stringr::str_split(gt.tidy$gt_AD, ",", simplify = TRUE)
gt.ad = data.table::as.data.table(gt.ad)

Avg.DP = base::mean(gt.tidy$gt_DP)
gt.ad$dp = gt.tidy$gt_DP

for(i in 1:ncol(gt.ad)) {
  
  gt.ad[[i]] = base::as.numeric(gt.ad[[i]])
  
}

gt.ad$pos = fix$POS

excluded.pos = gt.ad[which(gt.ad[[1]] == 0), ]$pos

ref = ref[which(!(ref$pos %in% as.numeric(excluded.pos))), ]

base::rm(vcf.file, reference.path)

fix = fix[base::which(fix$ALT != "<*>"), ]
fix$ALT = stringr::str_remove_all(fix$ALT, "\\,<\\*>")

ALT = stringr::str_split(fix$ALT, ",", simplify = TRUE)

out = base::list()

for(i in 1:base::ncol(ALT)) {
  
  who = base::which(ALT[,i] != "")
  
  out[[i]] = data.table::data.table(pos = fix[who, ]$POS,
                                    ref = fix[who, ]$REF,
                                    alt = ALT[who, i],
                                    ad = gt.ad[who, ][[i + 1]],
                                    dp = gt.ad[who, ]$dp)
  
}

fix = data.table::rbindlist(out)

base::rm(out, ALT, i, who)

comp = base::list()

exclude = c()

for(i in 1:base::nrow(fix)){
  
  ref.len = stringr::str_length(fix[i,]$ref)
  alt.len = stringr::str_length(fix[i,]$alt)
  
  if((ref.len == 1) & (alt.len == 1)) {
    
    temp = data.table::data.table(pos = fix[i,]$pos,
                                  logic.expr = "==",
                                  base = fix[i,]$alt)
    
    temp$whole.expr = base::paste(temp$pos, 
                                  temp$logic.expr, 
                                  "'", temp$base, "'", 
                                  sep = "")
    
    temp$ad = fix[i,]$ad
    temp$dp = fix[i,]$dp
    
    comp[[i]] = temp
    
  } else if((ref.len != 1) & (alt.len == 1)) {
    
    ref.split = stringr::str_split(fix[i,]$ref, "", simplify = TRUE)
    
    temp = data.table::data.table(pos = c(base::as.numeric(fix[i,]$pos), base::as.numeric(fix[i,]$pos) + 1:(ref.len - 1)),
                                  logic.expr = "==",
                                  base = c(fix[i,]$alt, base::rep("-", ref.len - 1)))
    
    temp$whole.expr = base::paste(temp$pos, 
                                  temp$logic.expr, 
                                  "'", temp$base, "'", 
                                  sep = "")
    
    temp$ad = fix[i,]$ad
    temp$dp = fix[i,]$dp
    
    comp[[i]] = temp
    
  } else if((ref.len != 1) & (alt.len != 1) & (ref.len == alt.len)) {
    
    ref.split = stringr::str_split(fix[i,]$ref, "", simplify = TRUE)
    alt.split = stringr::str_split(fix[i,]$alt, "", simplify = TRUE)
    
    ref.split = base::as.vector(ref.split)
    alt.split = base::as.vector(alt.split)
    
    diff = stringdist::stringdist(ref.split, alt.split)
    
    who = base::which(diff != 0)
    
    temp = data.table::data.table(pos = c(base::as.numeric(fix[i,]$pos) + who - 1),
                                  logic.expr = "==",
                                  base = alt.split[who])
    
    temp$whole.expr = base::paste(temp$pos, 
                                  temp$logic.expr, 
                                  "'", temp$base, "'", 
                                  sep = "")
    
    temp$ad = fix[i,]$ad
    temp$dp = fix[i,]$dp
    
    comp[[i]] = temp
    
  } else if((ref.len != 1) & (alt.len != 1) & (ref.len > alt.len)) {
    
    ref.split = stringr::str_split(fix[i,]$ref, "", simplify = TRUE)
    alt.split = stringr::str_split(fix[i,]$alt, "", simplify = TRUE)
    
    ref.split = base::as.vector(ref.split)
    alt.split = base::as.vector(alt.split)
    
    alt.j = 1
    
    gaps = c()
    
    for(j in 1:base::length(ref.split)) {
      
      if(ref.split[j] != alt.split[alt.j]) {
        
        gaps = base::c(gaps, 
                       base::as.numeric(fix[i,]$pos) + j - 1)
        
      } else {
        
        alt.j = base::ifelse((alt.j + 1) > alt.len, alt.len, alt.j + 1)
        
      }
      
    }
    
    if(base::is.null(gaps)) {
      
      gaps = base::as.numeric(fix[i,]$pos) + 1
      
    }
    
    temp = data.table::data.table(pos = gaps,
                                  logic.expr = "==",
                                  base = "-")
    
    temp$whole.expr = base::paste(temp$pos, 
                                temp$logic.expr, 
                                "'", temp$base, "'", 
                                sep = "")
    
    temp$ad = fix[i,]$ad
    temp$dp = fix[i,]$dp
    
    comp[[i]] = temp
    
    
  } else if((ref.len != 1) & (alt.len != 1) & (ref.len < alt.len)) {
    
    ref.split = stringr::str_split(fix[i,]$ref, "", simplify = TRUE)
    alt.split = stringr::str_split(fix[i,]$alt, "", simplify = TRUE)
    
    ref.split = base::as.vector(ref.split)
    alt.split = base::as.vector(alt.split)
    
    ref.j = ref.len
    
    insertions = base::c()
    
    for(j in base::length(alt.split):1) {
      
      if(alt.split[j] != ref.split[ref.j]) {
        
        insertions = base::c(insertions, j)
        
      } else {
        
        ref.j = base::ifelse((ref.j - 1) < 0, 1, ref.j - 1)
        
      }
      
    }
    
    if(length(insertions) == 1) {
      
      temp = data.table::data.table(pos = base::c(as.numeric(fix[i,]$pos) + insertions - 1),
                                    logic.expr = "==",
                                    base = alt.split[insertions])
      
      temp$whole.expr = base::paste(temp$pos, 
                                    temp$logic.expr, 
                                    "'", temp$base, "'", 
                                    sep = "")
      
      temp$ad = fix[i,]$ad
      temp$dp = fix[i,]$dp
      
      comp[[i]] = temp
      
    } else {
      
      exclude = base::c(exclude, i)
    }
  
  } else {
    
    exclude = base::c(exclude, i)
    
  }
  
}

comp = data.table::rbindlist(comp)

fix = comp

base::rm(list = base::setdiff(base::ls(), 
                              base::c("fix", "ref", "bases", "start.time", "decision.rules.path", "vcf.path", "nreads", "Avg.DP")))

fix$pos = as.numeric(fix$pos)

who = match(ref$pos, fix$pos)

ref$ad = fix[who, ]$dp - fix[who, ]$ad
ref$dp = fix[who, ]$dp 

ref[which(is.na(ref$dp)), ]$ad = 0 
ref[which(is.na(ref$dp)), ]$dp = 0 

ref$source = "ref"
fix$source = "var"

fix = base::rbind(ref, fix)

fix$pos = base::as.numeric(fix$pos)

fix = fix[base::order(fix$pos), ]


one.run = function(x, bases) {

  out = data.table::data.table(pos = base::as.numeric(x[1]),
                               logic.expr = "!=",
                               base = bases[!(bases %in% x[3])])
  
  out$whole.expr = base::paste(out$pos, 
                               out$logic.expr, 
                               "'", out$base, "'", 
                               sep = "")
  
  out$ad = x[5]
  out$dp = x[6]
  out$source = x[7]
  
  base::return(out)
}


comp = base::apply(fix, 1, one.run, bases = bases)

comp = data.table::rbindlist(comp)

comp = base::rbind(fix, comp)

comp = unique(comp)
comp = comp[order(comp$pos), ]

comp$ad = as.numeric(comp$ad)
comp$dp = as.numeric(comp$dp)

base::rm(bases, one.run, fix, ref)

# Create decision making rules ----------------------------------

decisions = data.table::fread(decision.rules.path, header = FALSE, sep = "\t", skip = 1)

rules = stringr::str_split(decisions$V2, ",", simplify = TRUE)
rules = data.table::as.data.table(rules)

# original.rules = rules
# 
# one.run = function(x, expr, AD, DP) {
#   
#   x = x[base::which(x != "")]
#   AD = AD[base::which(expr %in% x)]
#   AD = AD[which(AD != 0)]
#   
#   # DP = DP[base::which(expr %in% x)]
#   # DP = DP[which(DP != 0)]
#   
#   expr = unique(expr[base::which(expr %in% x)])
#   
#   
#   out = data.table::data.table(total.len = base::length(x),
#                                total.overlap = base::length(expr) # ,
#                                # avg.ad = base::mean(AD) # , 
#                                # avg.dp = base::mean(DP)
#                                )
#   
#   base::return(out)
#   
# }
# 
# rules.stats = base::apply(rules, 1, one.run, comp$whole.expr, comp$ad, comp$dp)
# rules.stats = data.table::rbindlist(rules.stats)
# 
# rules$n.rules = rules.stats$total.len
# 
# rules$n.totalOverlap = rules.stats$total.overlap

# rules$avg.ad = rules.stats$avg.ad
# rules$avg.dp = rules.stats$avg.dp

out = base::list()
k = 1

for(i in 1:(base::ncol(rules))) {
  
  who = rules[[i]] == ""
  
  if(base::length(which(who)) > 0){
    
    temp = decisions[base::which(who), ]
    # temp$total = rules[base::which(who), ]$n.rules
    temp$tree.overlap = i - 1
    # temp$total.overlap = rules[base::which(who), ]$n.totalOverlap
    
    
    # temp$total.avg.ad = rules[base::which(who), ]$avg.ad
    # temp$avg.dp = rules[base::which(who), ]$avg.dp
    
    out[[k]] = temp
    
    k = k + 1
    
    decisions = decisions[base::which(!who), ]
    rules = rules[base::which(!who), ]
    
  }
  
  if(base::nrow(decisions) > 0) {
    who = rules[[i]] %in% comp$whole.expr
    
    temp = decisions[base::which(!who), ]
    
    if(nrow(temp) > 0) {
      
      # temp$total = rules[base::which(!who), ]$n.rules
      temp$tree.overlap = i - 1
      # temp$total.overlap = rules[base::which(!who), ]$n.totalOverlap
      
      # temp$avg.ad = rules[base::which(!who), ]$avg.ad
      # temp$avg.dp = rules[base::which(!who), ]$avg.dp
      
      out[[k]] = temp
      
      k = k + 1
      
    }
    
    decisions = decisions[base::which(who), ]
    rules = rules[base::which(who), ]
    
  } else {
    
    break
    
  }
  
}

out = data.table::rbindlist(out)
 
one.run = function(x, expr, AD, DP, source) {
  
  rules = stringr::str_split(x[2], ",", simplify = TRUE)
  rules = as.vector(rules)

  tree.who = base::which(expr %in% rules[0:as.integer(x[3])])
  total.who = base::which(expr %in% rules)
  
  total.who.var = base::which(expr[which(source == "var")] %in% rules)
  
  tree.AD = AD[tree.who]
  tree.AD = tree.AD[which(tree.AD != 0)]
  
  total.AD = AD[total.who]
  total.AD = total.AD[which(total.AD != 0)]
  
  total.DP = DP[total.who]
  total.DP = total.DP[which(total.DP != 0)]

  expr.var = unique(expr[total.who.var])
  expr = unique(expr[total.who])
  


  stats = data.table::data.table(total.overlap = base::length(expr),
                                 total.overlap.var = base::length(expr.var),
                                 total.len = base::length(rules),
                                 tree.avg.ad = base::mean(tree.AD),
                                 total.avg.ad = base::mean(total.AD),
                                 total.avg.dp = base::mean(total.DP),
                                 total.sum.ad = base::sum(total.AD),
                                 total.sum.dp = base::sum(total.DP))

  base::return(stats)

}

rules.stats = base::apply(out, 1, one.run, comp$whole.expr, comp$ad, comp$dp, comp$source)
rules.stats = data.table::rbindlist(rules.stats)

# Generate output ----------------------------------

out$total.overlap = rules.stats$total.overlap
out$total.overlap.var = rules.stats$total.overlap.var
out$total = rules.stats$total.len
out$tree.avg.ad = rules.stats$tree.avg.ad
out$total.avg.ad = rules.stats$total.avg.ad
out$total.avg.dp = rules.stats$total.avg.dp
out$total.sum.ad = rules.stats$total.sum.ad
out$total.sum.dp = rules.stats$total.sum.dp

out$tree.ratio = out$tree.overlap / out$total
out$total.ratio = out$total.overlap / out$total
out$total.ratio.var = out$total.overlap.var / out$total

out$avg.dp = Avg.DP
out$total.run.reads = nreads

out[which(is.na(out$tree.avg.ad)), ]$tree.avg.ad = 0
out[which(is.na(out$total.avg.ad)), ]$total.avg.ad = 0
out[which(is.na(out$total.avg.dp)), ]$total.avg.dp = 0
out[which(is.na(out$total.sum.ad)), ]$total.sum.ad = 0
out[which(is.na(out$total.sum.dp)), ]$total.sum.dp = 0

out = out[base::order(out$tree.ratio, out$total.ratio, decreasing = TRUE), ]

base::colnames(out) = base::c("Lineage", "Rules", 
                              "Tree.Overlap", "Total.Overlap", "Total.Overlap.Var", "Total",
                              "Tree.Avg.AD", "Total.Avg.AD", "Total.Avg.DP",
                              "Total.Sum.AD", "Total.Sum.DP",
                              "Tree.Ratio", "Total.Ratio", "Total.Ratio.Var",
                              "Avg.DP", "Total.Run.Reads")

out = out[,c("Lineage", 
             "Rules", 
             "Total", 
             "Tree.Overlap", 
             "Total.Overlap", 
             "Total.Overlap.Var",
             "Tree.Ratio", 
             "Total.Ratio", 
             "Total.Ratio.Var",
             "Tree.Avg.AD", 
             "Total.Avg.AD", 
             "Total.Avg.DP",
             "Total.Sum.AD", 
             "Total.Sum.DP",
             "Avg.DP", 
             "Total.Run.Reads")]

out$Sig = out$Total.Sum.AD / out$Total.Sum.DP

out[which(is.na(out$Sig)), ]$Sig = 0

end.time = base::Sys.time()

utils::write.table(out, 
                   paste(stringr::str_replace(vcf.path, ".vcf", "-lineagespot.tsv"), sep = ""), 
                   row.names = FALSE, 
                   quote = FALSE, 
                   sep = "\t")

base::rm(list = base::setdiff(base::ls(), 
                              base::c("out", "comp", "decisions", "start.time", "end.time")))














