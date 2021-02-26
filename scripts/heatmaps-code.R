# Libraries needed ----------------------------------

library(vcfR)
library(data.table)
library(stringr)
library(seqinr)
library(stringdist)
library(dplyr)
library(saveImageHigh)
library(stringi)
library(tidyr)


base::rm(list = ls())


# Inputs ----------------------------------
timepoints <- c("02-14/12/2020","05-11/02/2021", "12-18/02/2021")

bases = base::c("A", "G", "C", "T", "-")

reference.path = "ref/NC_045512.fasta"

files <- list.files(path = "./vcf-files", pattern = "freebayes")

decision.rules.path = "ref/decision_tree_rules_2021_02_17.txt"

utils::download.file(
  url = "https://raw.githubusercontent.com/cov-lineages/pangoLEARN/master/pangoLEARN/data/decision_tree_rules.txt",
  destfile = decision.rules.path,
  method = "curl"
)

# Create reference based rules ----------------------------------

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
list.heatmap <- list()
k <- 1

for (vcf in files){
  
  vcf.file = vcfR::read.vcfR(paste0("vcf-files/",vcf), verbose = FALSE)
  fix = data.table::as.data.table(vcf.file@fix)
  
  
  info <- fix[,c("POS","INFO")]
  info$INFO <- str_split(info$INFO,"ANN", simplify = TRUE)[,2]
  
  info <- info %>% separate(INFO,c("Allele", "Annotation",
                                    "Annotation_Impact", "Gene_Name" ,
                                     "Gene_ID", "Feature_Type",
                                     "Feature_ID ","Transcript_BioType","Rank",
                                     "HGVS.c" , "HGVS.p" ,"cDNA.pos_cDNA.length",
                                      "CDS.pos_CDS.length", "AA.pos_AA.length",
                                      "Distance", "ERRORS_WARNINGS_INFO" ),"\\|")
  info <- info[,c("POS","Gene_Name","HGVS.c","HGVS.p")]
  
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
  
  #base::rm(vcf.file, reference.path)
  
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
  
  fix <- full_join(fix,info, by = c("pos" ="POS"))
  
  #base::rm(out, ALT, i, who)
  
  comp = base::list()
  
  exclude = c()
  
  for(i in 1:base::nrow(fix)){
    
    ref.len = stringr::str_length(fix[i,]$ref)
    alt.len = stringr::str_length(fix[i,]$alt)
    
    if((ref.len == 1) & (alt.len == 1)) {
      
      temp = data.table::data.table(pos = fix[i,]$pos,
                                    rule = paste0(fix[i,]$ref," -> ",fix[i,]$alt),
                                    freq = fix[i,]$ad / fix[i,]$dp,
                                    Gene_Name = fix[i,]$Gene_Name,
                                    HGVS.c = fix[i,]$HGVS.c,
                                    HGVS.p = fix[i,]$HGVS.p)
      
      
      comp[[i]] = temp
      
    } else if((ref.len != 1) & (alt.len == 1)) {
      
      ref.split = stringr::str_split(fix[i,]$ref, "", simplify = TRUE)
      
      temp = data.table::data.table(pos = c(base::as.numeric(fix[i,]$pos), 
                                            base::as.numeric(fix[i,]$pos) + 1:(ref.len - 1)),
                                    rule = paste0(c(ref.split)[-1]," -> -"),
                                    freq = fix[i,]$ad / fix[i,]$dp,
                                    Gene_Name = fix[i,]$Gene_Name,
                                    HGVS.c = fix[i,]$HGVS.c,
                                    HGVS.p = fix[i,]$HGVS.p)
  
      
      comp[[i]] = temp
      
    } else if((ref.len != 1) & (alt.len != 1) & (ref.len == alt.len)) {
      
      ref.split = stringr::str_split(fix[i,]$ref, "", simplify = TRUE)
      alt.split = stringr::str_split(fix[i,]$alt, "", simplify = TRUE)
      
      ref.split = base::as.vector(ref.split)
      alt.split = base::as.vector(alt.split)
      
      diff = stringdist::stringdist(ref.split, alt.split)
      
      who = base::which(diff != 0)
      
      temp = data.table::data.table(pos = c(base::as.numeric(fix[i,]$pos) + who - 1),
                                    rule = paste0(c(ref.split)[who]," -> ",c(alt.split)[who]),
                                    freq = fix[i,]$ad / fix[i,]$dp,
                                    Gene_Name = fix[i,]$Gene_Name,
                                    HGVS.c = fix[i,]$HGVS.c,
                                    HGVS.p = fix[i,]$HGVS.p)
      
      
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
      
      index <- gaps - as.numeric(fix$pos[i]) + 1
      
      temp = data.table::data.table(pos = gaps,
                                    rule = paste0(ref.split[index]," -> -"),
                                    freq = fix[i,]$ad / fix[i,]$dp,
                                    Gene_Name = fix[i,]$Gene_Name,
                                    HGVS.c = fix[i,]$HGVS.c,
                                    HGVS.p = fix[i,]$HGVS.p)
      
   
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
                                      rule = paste0("- -> ",alt.split[insertions]),
                                      freq = fix[i,]$ad / fix[i,]$dp,
                                      Gene_Name = fix[i,]$Gene_Name,
                                      HGVS.c = fix[i,]$HGVS.c,
                                      HGVS.p = fix[i,]$HGVS.p)
       
        
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
  #fix = group_by(fix,pos,rule) %>%
  #  summarise(pos,rule, freq = mean(freq)) %>%
  #  as.data.table() %>%
  #  unique()
  
  #fix <- group_by(fix,pos) %>%
  #  summarise(pos,rule = paste(rule,collapse =" | "), freq = mean(freq)) %>%
  #  as.data.table() %>%
  #  unique()
  
  fix$ID <- paste0(fix$pos, "  ", fix$rule)
  fix <- fix[,c("ID","pos","freq", "Gene_Name", "HGVS.c"  ,  "HGVS.p" )]
  
  fix = group_by(fix,ID) %>%
    summarise(ID, freq = mean(freq),pos,Gene_Name ,HGVS.c, HGVS.p) %>%
    as.data.table() %>%
    unique()
  
  fix$pos <- as.numeric(fix$pos)
  
  colnames(fix) <- c("ID",paste0(colnames(fix)[-1],"_",files[k]))
  list.heatmap[[k]] <- fix
  k <- k + 1
  
}

names(list.heatmap) <- timepoints

common_join <- list.heatmap[[2]]

for (l in c(1,3)){
  
  temp <- list.heatmap[[l]]
  common_join <- inner_join(common_join,temp[,1:2])
}

common_join$pos_L1_S22_L001_freebayes.ann.vcf <- as.numeric(common_join$pos_L1_S22_L001_freebayes.ann.vcf)

common_join <- common_join[order(common_join[,9], decreasing = F),]


data.common <- common_join[ ,c(2,7,8)]
data.common <- as.matrix(data.common)
row.names(data.common) <- as.character(common_join$`pos_Sewage-L2_S10_L001_freebayes.ann.vcf`)
colnames(data.common) <- timepoints[c(2,1,3)]
#data.common <- data.common[order(row.names(data.common), decreasing = F),]

#columns <- c(seq(2,2*length(files),2))
#rules.common <- common_join[,..columns]
#row.names(rules.common) <- as.character(common_join$pos)
#rules.common <- rules.common[order(row.names(rules.common), decreasing = T),]

library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c( 0,1), c("#FFF6FF", "red"))
#col_fun2 = colorRamp2( unique(annot),c("red", "green","blue", "yellow" ,"grey","pink","brown"))
#col_fun(seq(-2, 2))

annot <- str_split(common_join$ID, pattern ="  ", n=2)
annot <- unlist(annot)
annot <- annot[c(F,T)]

genes <- common_join$`Gene_Name_Sewage-L2_S10_L001_freebayes.ann.vcf`

HGVS.c <- common_join$`HGVS.c_Sewage-L2_S10_L001_freebayes.ann.vcf`

HGVS.p <- common_join$`HGVS.p_Sewage-L2_S10_L001_freebayes.ann.vcf`

extra = rowAnnotation(variants = anno_text(annot, location = 0.5, just = "center",
                                              gp = gpar(fontsize = 11, fill = "#DCE5F3", col = "black", border = "white"),
                                              width = max_text_width(annot)*1.2),
                      Gene_name = anno_text(genes,
                                            location = 0.5, just = "center",
                                            gp = gpar(fontsize = 11, fill = "#F5CDF7", col = "black", border = "white"),
                                            width = max_text_width(genes)*1.2),
                      HGVS.c = anno_text(HGVS.c,
                                            location = 0.5, just = "center",
                                            gp = gpar(fontsize = 11, fill = "#D3F7CD", col = "black", border = "white"),
                                            width = max_text_width(HGVS.c)),
                      HGVS.p = anno_text(HGVS.p,
                                         location = 0.5, just = "center",
                                         gp = gpar(fontsize = 11, fill = "#F8E5A6", col = "black", border = "white"),
                                         width = max_text_width(HGVS.p))
                      ) 
                      #annotation_name_side = "top", annotation_name_rot = 0)
#"T -> C" "C -> T" "C -> A" "T -> G" "G -> T" "A -> G" "A -> C"

#row_ha = rowAnnotation(variants = anno_simple(annot,width = unit(3, "cm")))

data.common <- data.common[,c(2,1,3)]

graph1 <-Heatmap(data.common, name = "frequency", 
                  col = col_fun,
                  heatmap_width = unit(21, "cm"), heatmap_height = unit(8, "cm"),
                  border = T,
                  #border_gp = gpar(col = "grey"),
                  rect_gp = gpar(col = "white", lwd = 1),
                  cluster_rows = FALSE,cluster_columns  = FALSE,
                 column_names_gp = gpar(fontsize = 13),
                 row_names_gp = gpar(fontsize = 14),
                  column_names_rot = 50,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.2f", data.common[i,j]), x, y, gp = gpar(fontsize = 12))
                  },
                 heatmap_legend_param = list(
                   labels_gp = gpar(fontsize = 13),
                   title_gp =gpar(fontsize = 15),
                   legend_height = unit(2.5, "cm"),
                   legend_width = unit(4, "cm")
                 ),
                 right_annotation = extra,
                 row_split = factor(genes, levels = unique(genes)), #factor(rep(LETTERS[1:3], 6), levels = LETTERS[3:1])
                 row_order = order(as.numeric(row.names(data.common)), decreasing = F),row_title_rot = 0,
                 cluster_row_slices = F,
                 )

#print(graph1)

save_image(print(graph1),"common variants.png", width = 10, height = 10)


threshold <- 0.3

data.threshold <- list.heatmap[[1]]
data.threshold <- data.threshold[which(data.threshold[,2] >= threshold), ]


for (l in 2:length(files)){
  temp <- list.heatmap[[l]]
  temp <- temp[which(temp[,2] >= threshold), ]
  data.threshold <- full_join(data.threshold,temp)
}

data.threshold[is.na(data.threshold)] <- ""
data.threshold <- data.threshold %>%
                  group_by(ID) %>%
                  mutate(genes = max(Gene_Name_L1_S22_L001_freebayes.ann.vcf,
                                     `Gene_Name_Sewage-L2_S10_L001_freebayes.ann.vcf`,
                                     Gene_Name_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf, 
                                     na.rm=TRUE))

data.threshold <- data.threshold %>%
  group_by(ID) %>%
  mutate(HGVS.p = max(HGVS.p_L1_S22_L001_freebayes.ann.vcf,
                      `HGVS.p_Sewage-L2_S10_L001_freebayes.ann.vcf`,
                      HGVS.p_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf, 
                     na.rm=TRUE))
 
positions <- str_split(data.threshold$ID, pattern ="  ", n=2)
positions <- unlist(positions)
positions <- positions[c(T,F)]

data.threshold$positions <- as.numeric(positions)

data.threshold <- data.threshold[order(data.threshold$positions, decreasing = F),]

columns <- c(2,7,12)
data.graph2 <- data.threshold[ ,c(2,7,12)]
data.graph2 <- as.matrix(data.graph2)


#data.graph2 <- data.graph2[order(colnames(data.graph2), decreasing = T),]

data.graph2 <- apply(data.graph2,2,function(x) as.numeric(x))
data.graph2[is.na(data.graph2)] <- 0
colnames(data.graph2) <- timepoints 
row.names(data.graph2) <- data.threshold$positions

extra2 = rowAnnotation(Gene_name = anno_text(data.threshold$genes,
                                            location = 0.5, just = "center",
                                            gp = gpar(fontsize = 6, fill = "#F5CDF7", col = "black", border = "white"),
                                            width = unit(1.5, "cm")),
                      HGVS.p = anno_text(data.threshold$HGVS.p,
                                         location = 0.5, just = "center",
                                         gp = gpar(fontsize = 6, fill = "#F8E5A6", col = "black", border = "white"),
                                         width = unit(3, "cm"))
) 




graph2 <-Heatmap(data.graph2, name = "frequency", 
                 col = col_fun,
                 heatmap_width = unit(11, "cm"), heatmap_height = unit(1200, "cm"),
                 border = T,
                 #border_gp = gpar(col = "grey"),
                 rect_gp = gpar(col = "white", lwd = 1),
                 cluster_rows = FALSE,cluster_columns  = FALSE,
                 row_order = order(as.numeric(row.names(data.graph2)), decreasing = F),
                 column_names_gp = gpar(fontsize = 10),
                 row_names_gp = gpar(fontsize = 7),
                 column_names_rot = 50,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", data.graph2[i,j]), x, y, gp = gpar(fontsize = 6))
                 },
                 heatmap_legend_param = list(
                    labels_gp = gpar(fontsize = 9),
                    title_gp =gpar(fontsize = 10),
                   legend_height = unit(3, "cm"),
                   legend_width = unit(4, "cm"),
                   at = c(0, threshold,1)
                 ),
                 right_annotation = extra2,
                 row_split = factor(data.threshold$genes, levels = unique(data.threshold$genes)), 
                 row_title_gp = gpar(fontsize = 9), row_title_rot = 0,
                 cluster_row_slices = F)

save_image(print(graph2),"most frequent variants.png", width = 10, height = 15)

############### All variants #############
#data.threshold <- list.heatmap[[1]]

#for (l in 2:length(files)){
#  temp <- list.heatmap[[l]]
#  data.threshold <- full_join(data.threshold,temp)
#}

#data.threshold[is.na(data.threshold)] <- ""
#data.threshold <- data.threshold %>%
#  group_by(ID) %>%
#   mutate(genes = max(Gene_Name_L1_S22_L001_freebayes.ann.vcf,
#                      `Gene_Name_Sewage-L2_S10_L001_freebayes.ann.vcf`,
#                      Gene_Name_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf, 
#                      na.rm=TRUE))
# 
# data.threshold <- data.threshold %>%
#   group_by(ID) %>%
#   mutate(HGVS.p = max(HGVS.p_L1_S22_L001_freebayes.ann.vcf,
#                       `HGVS.p_Sewage-L2_S10_L001_freebayes.ann.vcf`,
#                       HGVS.p_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf, 
#                       na.rm=TRUE))
# 
# data.threshold <- data.threshold[,c(1,2,7,12,17,18)]
# for (one.gene in unique(data.threshold$genes)){
#   
#   one.table <- data.threshold[which(data.threshold$genes == one.gene),]
#   
#   positions <- str_split(one.table$ID, pattern ="  ", n=2)
#   positions <- unlist(positions)
#   positions <- positions[c(T,F)]
#   
#   one.table$positions <- as.numeric(positions)
#   
#   one.table<- one.table[order(one.table$positions, decreasing = F),]
# 
#   data.graph2 <- one.table[ ,c(2:4)]
#   data.graph2 <- as.matrix(data.graph2)
#   
#   
#   #data.graph2 <- data.graph2[order(colnames(data.graph2), decreasing = T),]
#   
#   data.graph2 <- apply(data.graph2,2,function(x) as.numeric(x))
#   data.graph2[is.na(data.graph2)] <- 0
#   colnames(data.graph2) <- timepoints 
#   row.names(data.graph2) <- one.table$positions
#   
#   extra2 = rowAnnotation(Gene_name = anno_text(one.table$genes,
#                                                location = 0.5, just = "center",
#                                                gp = gpar(fontsize = 6, fill = "#F5CDF7", col = "black", border = "white"),
#                                                width = unit(1.5, "cm")),
#                          HGVS.p = anno_text(one.table$HGVS.p,
#                                             location = 0.5, just = "center",
#                                             gp = gpar(fontsize = 6, fill = "#F8E5A6", col = "black", border = "white"),
#                                             width = unit(3, "cm"))
#   ) 
#   
#   
#   
#   
#   graph2 <-Heatmap(data.graph2, name = "frequency", 
#                    col = col_fun,
#                    heatmap_width = unit(11, "cm"), heatmap_height = unit(900, "cm"),
#                    border = T,
#                    #border_gp = gpar(col = "grey"),
#                    rect_gp = gpar(col = "white", lwd = 1),
#                    cluster_rows = FALSE,cluster_columns  = FALSE,
#                    row_order = order(as.numeric(row.names(data.graph2)), decreasing = F),
#                    column_names_gp = gpar(fontsize = 10),
#                    row_names_gp = gpar(fontsize = 7),
#                    column_names_rot = 50,
#                    cell_fun = function(j, i, x, y, width, height, fill) {
#                      grid.text(sprintf("%.2f", data.graph2[i,j]), x, y, gp = gpar(fontsize = 6))
#                    },
#                    heatmap_legend_param = list(
#                      labels_gp = gpar(fontsize = 9),
#                      title_gp =gpar(fontsize = 10),
#                      legend_height = unit(3, "cm"),
#                      legend_width = unit(4, "cm"),
#                      at = c(0, 1)
#                    ),
#                    right_annotation = extra2,
#                    #row_split = factor(data.threshold$genes, levels = unique(data.threshold$genes)), 
#                    #row_title_gp = gpar(fontsize = 9), row_title_rot = 0,
#                    cluster_row_slices = F)
#   
#   #saveImageHigh::save_image(print(graph2),paste0(one.gene,"_variants.png"), width = 5, height = 400, res = 380)
#   save_as_pdf(print(graph2),paste0(one.gene,"_variants.pdf"), width = 6, height = 900)
# }

########### Lineages heatmaps ############
#lineages <- read.csv("Sewage-L2_S10_L001_freebayes.ann-lineagespot-overrules-v2.tsv", sep ="\t")
#l.rules <- lineages[,c("Lineage","Overlap.Rules", "Rules")]
#l.rules <- l.rules[which(l.rules$Lineage %in% c("B.1.1.7","B.1.351","B.1.177")),]

uk <- c(6954,5388, 3267, 28977, 28111, 28048, 27972, 24914, 24506, 23709, 23604, 23271, 23063, 
                                28280, 28281, 28282, 21991,21992,21993, 21765,21766, 21767, 21768,21769, 21770,
                                 11288, 11289, 11290, 11291, 11292, 11293, 11294, 11295, 11296)

africa <- c(28887, 26156, 25904, 25563, 23664, 23063, 23012, 22813, 21801, 10323, 5230, 1059)

reference.path = "ref/NC_045512.fasta"

ref = seqinr::read.fasta(reference.path)
ref = base::as.character(ref[[1]])

ref = data.table::data.table(pos = 1:base::length(ref),
                             logic.expr = "==",
                             base = stringr::str_to_upper(ref))
k <- 1
uk.list <- list()
africa.list <- list()

for (vcf in files){
  
  #print(vcf)
  vcf.file = vcfR::read.vcfR(paste0("vcf-files/",vcf), verbose = FALSE)
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
  
  #base::rm(vcf.file, reference.path)
  
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
  fix$pos <- as.numeric(fix$pos)
  
  uk.fix <- fix[which(fix$pos %in% uk),]
  who_23604 <- which((uk.fix$pos == 23604) & str_sub(uk.fix$ref,1,1) == str_sub(uk.fix$alt,1,1))
  if (length(who_23604) != 0){
    uk.fix <- uk.fix[-c(who_23604),]
  }
  
  who_23604 <- which((uk.fix$pos == 23604))
  uk.fix$ref[who_23604] <- str_sub(uk.fix$ref[who_23604],1,1)
  uk.fix$alt[who_23604] <- str_sub(uk.fix$alt[who_23604],1,1)
  
  who_24506 <- which((uk.fix$pos == 24506))
  uk.fix$ref[who_24506] <- str_sub(uk.fix$ref[who_24506],1,1)
  uk.fix$alt[who_24506] <- str_sub(uk.fix$alt[who_24506],1,1)
    
  uk.fix <- group_by(uk.fix,pos) %>%
    summarise(pos,ref,alt = paste(alt, collapse = ","), ad = sum(ad), dp) %>%
    as.data.table() %>%
    unique()
  
  
  uk.list[[k]] <- uk.fix
  
  africa.fix <- fix[which(fix$pos %in% africa),]
  
  africa.fix <- group_by(africa.fix,pos) %>%
    summarise(pos,ref,alt = paste(alt, collapse = ","), ad = sum(ad), dp) %>%
    as.data.table() %>%
    unique()
  africa.list[[k]] <- africa.fix
  k <- k+1
}

for (l in 1:length(uk.list)){
  temp <- data.table(pos = uk,af = 0, rules = "")
  temp$pos <- temp$pos
  one.table <- uk.list[[l]]
  
  temp <- full_join(temp,one.table)
  
  temp$af<- temp$ad/temp$dp
  temp$rules<- paste0(temp$ref," -> ",temp$alt)
  
  temp <- temp[,c("pos","af","rules")]
  
  vcf.file = vcfR::read.vcfR(paste0("vcf-files/",files[l]), verbose = FALSE)
  fix = data.table::as.data.table(vcf.file@fix)
  
  who.present <- which((temp$pos %in% fix$POS) & temp$af != NA)
  
  temp$af[who.present] <- 0
  temp$rules <- str_replace(temp$rules,"NA -> NA","")
  colnames(temp) <- c("pos",paste0("af_",files[l]),paste0("rules_",files[l]))
  
  uk.list[[l]] <- temp 
}

for (l in 1:length(africa.list)){
  temp <- data.table(pos = africa,af = 0, rules = "")
  temp$pos <- temp$pos
  one.table <- africa.list[[l]]
  
  temp <- full_join(temp,one.table)
  
  temp$af<- temp$ad/temp$dp
  temp$rules<- paste0(temp$ref," -> ",temp$alt)
  
  temp <- temp[,c("pos","af","rules")]
  
  vcf.file = vcfR::read.vcfR(paste0("vcf-files/",files[l]), verbose = FALSE)
  fix = data.table::as.data.table(vcf.file@fix)
  
  who.present <- which((temp$pos %in% fix$POS) & temp$af != NA)
  
  temp$af[who.present] <- 0
  temp$rules <- str_replace(temp$rules,"NA -> NA","")
  
  colnames(temp) <- c("pos",paste0("af_",files[l]),paste0("rules_",files[l]))
  africa.list[[l]] <- temp
}

merged.uk <- inner_join(uk.list[[1]],uk.list[[2]])
merged.uk <- inner_join(merged.uk,uk.list[[3]])

merged.africa <- inner_join(africa.list[[1]],africa.list[[2]])
merged.africa <- inner_join(merged.africa,africa.list[[3]])

who_28280 <- which(merged.uk$pos %in% seq(28280,28282,1))
merged.uk$pos[who_28280] <- "[28280-28282]"

who_21991 <- which(merged.uk$pos %in% seq(21991,21993,1))
merged.uk$pos[who_21991] <- "[21991-21993]"

who_21765 <- which(merged.uk$pos %in% seq(21765,21770,1))
merged.uk$pos[who_21765] <- "[21765-21770]"

who_11288 <- which(merged.uk$pos %in% seq(11288,11296,1))
merged.uk$pos[who_11288] <- "[11288-11296]"

merged.uk[is.na(merged.uk)] <- 0

merged.uk <- group_by(merged.uk,pos) %>%
  summarise(pos,af_L1_S22_L001_freebayes.ann.vcf = sum(af_L1_S22_L001_freebayes.ann.vcf),
            rules_L1_S22_L001_freebayes.ann.vcf = paste(rules_L1_S22_L001_freebayes.ann.vcf, collapse = " "),
            `af_Sewage-L2_S10_L001_freebayes.ann.vcf` = sum(`af_Sewage-L2_S10_L001_freebayes.ann.vcf`),
            `rules_Sewage-L2_S10_L001_freebayes.ann.vcf` = paste(`rules_Sewage-L2_S10_L001_freebayes.ann.vcf`, collapse= " "),
            af_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf = sum(af_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf),
            rules_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf = paste(rules_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf, collapse = " "),) %>%
  as.data.table() %>%
  unique()

columns <- c(2,4,6)
data.graph3 <- merged.uk[ ,..columns]
data.graph3 <- as.matrix(data.graph3)

data.graph3 <- apply(data.graph3,2,function(x) as.numeric(x))
#data.graph2[is.na(data.graph2)] <- 0
colnames(data.graph3) <- timepoints 
row.names(data.graph3) <- merged.uk$pos

extra3 = rowAnnotation(variants_1 = anno_text(merged.uk$rules_L1_S22_L001_freebayes.ann.vcf,
                                             location = 0.5, just = "center",
                                             gp = gpar(fontsize = 6, fill = "#DCE5F3", col = "black", border = "white"),
                                             width = unit(1, "cm")),
                       variants_2 = anno_text(merged.uk$`rules_Sewage-L2_S10_L001_freebayes.ann.vcf`,
                                          location = 0.5, just = "center",
                                          gp = gpar(fontsize = 6, fill = "#DCE5F3", col = "black", border = "white"),
                                          width = unit(1.5, "cm")),
                       variants_3 = anno_text(merged.uk$rules_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf,
                                              location = 0.5, just = "center",
                                              gp = gpar(fontsize = 6, fill = "#DCE5F3", col = "black", border = "white"),
                                              width = unit(1.5, "cm"))
) 


data.graph3[is.na(data.graph3)] <- 0

graph3 <-Heatmap(data.graph3, name = "frequency", 
                 col = col_fun,
                 heatmap_width = unit(10, "cm"), heatmap_height = unit(9, "cm"),
                 border = T,
                 #border_gp = gpar(col = "grey"),
                 rect_gp = gpar(col = "white", lwd = 1),
                 cluster_rows = FALSE,cluster_columns  = FALSE,
                 row_order = order(row.names(data.graph3), decreasing = T),
                 column_names_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 7),
                 column_names_rot = 50,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", data.graph3[i,j]), x, y, gp = gpar(fontsize = 6.5))
                 },
                 heatmap_legend_param = list(
                   labels_gp = gpar(fontsize = 7),
                   title_gp =gpar(fontsize = 8),
                   legend_height = unit(3, "cm"),
                   legend_width = unit(1, "cm")
                 ),
                 right_annotation = extra3,
                 #row_split = factor(data.threshold$genes, levels = unique(data.threshold$genes)), 
                 #row_title_gp = gpar(fontsize = 9), row_title_rot = 0,
                 #cluster_row_slices = F
                 )

save_image(print(graph3),"uk heatmap.png", width = 10, height = 10)



columns <- c(2,4,6)
data.graph4 <- merged.africa[ ,..columns]
data.graph4 <- as.matrix(data.graph4)

data.graph4 <- apply(data.graph4,2,function(x) as.numeric(x))
#data.graph2[is.na(data.graph2)] <- 0
colnames(data.graph4) <- timepoints 
row.names(data.graph4) <- merged.africa$pos

extra4 = rowAnnotation(variants_1 = anno_text(merged.africa$rules_L1_S22_L001_freebayes.ann.vcf,
                                              location = 0.5, just = "center",
                                              gp = gpar(fontsize = 6, fill = "#DCE5F3", col = "black", border = "white"),
                                              width = unit(1, "cm")),
                       variants_2 = anno_text(merged.africa$`rules_Sewage-L2_S10_L001_freebayes.ann.vcf`,
                                              location = 0.5, just = "center",
                                              gp = gpar(fontsize = 6, fill = "#DCE5F3", col = "black", border = "white"),
                                              width = unit(1, "cm")),
                       variants_3 = anno_text(merged.africa$rules_Sewage_CoV19_L3_S1_L001_freebayes.ann.vcf,
                                              location = 0.5, just = "center",
                                              gp = gpar(fontsize = 6, fill = "#DCE5F3", col = "black", border = "white"),
                                              width = unit(1, "cm"))
) 


data.graph4[is.na(data.graph4)] <- 0

graph4 <-Heatmap(data.graph4, name = "frequency", 
                 col = col_fun,
                 heatmap_width = unit(8, "cm"), heatmap_height = unit(8, "cm"),
                 border = T,
                 #border_gp = gpar(col = "grey"),
                 rect_gp = gpar(col = "white", lwd = 1),
                 cluster_rows = FALSE,cluster_columns  = FALSE,
                 row_order = order(as.numeric(row.names(data.graph4)), decreasing = F),
                 column_names_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 7),
                 column_names_rot = 50,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f", data.graph4[i,j]), x, y, gp = gpar(fontsize = 6.5))
                 },
                 heatmap_legend_param = list(
                   labels_gp = gpar(fontsize = 7),
                   title_gp =gpar(fontsize = 8),
                   legend_height = unit(3, "cm"),
                   legend_width = unit(1, "cm")
                 ),
                 right_annotation = extra4,
                 #row_split = factor(data.threshold$genes, levels = unique(data.threshold$genes)), 
                 #row_title_gp = gpar(fontsize = 9), row_title_rot = 0,
                 #cluster_row_slices = F
)

save_image(print(graph4),"africa heatmap.png", width = 10, height = 10)
