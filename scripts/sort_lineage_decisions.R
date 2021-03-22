library(data.table)
library(stringr)

rm(list = ls())

decisions = fread("decision_tree_rules_2021_01_13.txt", header = FALSE)

# comp = as.data.table(str_split(decisions$V2, ",", simplify = TRUE))
# 
# comp.pos = matrix(data = 0, nrow = nrow(comp), ncol = ncol(comp))
# comp.pos = as.data.table(comp.pos)
# 
# 
# for(i in 1:ncol(comp)) {
#   
#   
#   
#   temp = str_split(comp[[i]], "='", simplify = TRUE)[,1]
#   temp = str_sub(temp, 1, -2)
#   
#   comp.pos[[i]] = as.numeric(temp)
#   
#   cat(c(i, unique(temp), "\n"))
#   
# }
# 
# rm(temp, i)

one.run = function(x) {
  
  comp = str_split(x[2], ",", simplify = TRUE)
  out = str_split(comp, "='", simplify = TRUE)[,1]
  out = str_sub(out, 1, -2)
  out = as.numeric(out)

  comp = comp[order(out)]

  return(data.table(lineage = x[1],
                    rules = paste(comp, collapse = ",")))
  
  
}

decisions.sort = apply(decisions, 1, one.run)

decisions.sort = rbindlist(decisions.sort)
decisions.sort = decisions.sort[order(decisions.sort$lineage), ]

write.table(decisions.sort, file = "decision_tree_rules_2021_01_13_sorted.txt", sep = "\t", row.names = FALSE, quote = FALSE)

