check_deletions <- function(input_data,
                            deletions,
                            timepoints,
                            aa_abbreviations){
  
  input_data$POS <- as.numeric(input_data$POS)
  input_data$match <- input_data$HGVS.p
  
  # create deletion string in the same format
  for (i in c(1:nrow(aa_abbreviations))){
    
    input_data$match <- str_remove_all(input_data$match, aa_abbreviations$Three_Letter[i])
  }
  
  input_data$match <- str_remove_all(input_data$match, "p.")
  
  deletions$match <- paste0(deletions$codon_num, "_", deletions$codon_end, "del")
  who_one <- which(deletions$codon_num == deletions$codon_end)
  deletions$match[who_one] <- paste0(as.numeric(deletions$codon_num[who_one]) + 1, "del")
  
  # find matches
  specific_data <- input_data[which(input_data$match %in% deletions$match), ]
  specific_data <- specific_data[which(specific_data$TYPE == "del"), ]
  
  specific_data <- specific_data[order(specific_data$POS, decreasing = F), ]
  
  data_specific <- matrix(0, nrow = nrow(deletions), ncol = nrow(timepoints))
  row.names(data_specific) <- deletions$match
  colnames(data_specific) <- timepoints$samples
  
  #initialize annotations
  annotations <- data.table(annot = character(nrow(deletions)),
                            genes = character(nrow(deletions)),
                            HGVS.c = character(nrow(deletions)),
                            HGVS.p = character(nrow(deletions)),
                            POS = character(nrow(deletions)),
                            match = character(nrow(deletions)))
  
  row.names(annotations) <- deletions$match
  
  for (k in timepoints$samples){
    
    temp <- specific_data[which(specific_data$sample == k), 
                          c("POS", "REF", "ALT", "TYPE", "Gene_Name", "HGVS.c", "HGVS.p", "AF", "match")]
    
    # who.one <- which(deletions$change_length_nt == 1)
    # one.replacement <- deletions$POS[who.one]
    # not.really.ones <- which((temp$POS %in% one.replacement) & 
    #                            ((nchar(temp$REF))-(nchar(temp$ALT)) != 1))
    # 
    # who.many <- which(deletions$change_length_nt > 1)
    # many.replacement <- deletions$POS[who.many]
    # not.really.manys <- which((temp$POS %in% many.replacement) & 
    #                             ((nchar(temp$REF))-(nchar(temp$ALT)) != deletions$change_length_nt))
    # 
    # removed <-c(not.really.ones, not.really.manys)
    # 
    # if (length(removed > 0)){
    #   temp <- temp[-removed, ]
    # }
    
    
    temp <- group_by(temp, POS) %>%
      summarise(POS, REF,
                ALT = paste(ALT, collapse = ","), Gene_Name, HGVS.c, HGVS.p ,
                AF = sum(AF), match) %>%
      as.data.table %>%
      unique()
    
    
    data_specific[as.character(temp$match), k] <- temp$AF
    
    to.be.annotated <- which(row.names(annotations) %in% as.character(temp$match))
    
    annotations[to.be.annotated, c("annot")] <- paste0(temp$REF, " -> ", temp$ALT)
    annotations$genes[to.be.annotated] <- temp$Gene_Name
    annotations$HGVS.c[to.be.annotated] <- temp$HGVS.c
    annotations$HGVS.p[to.be.annotated] <- temp$HGVS.p
    annotations$POS[to.be.annotated] <- temp$POS
    annotations$match[to.be.annotated] <- temp$match
    
    who.huge <- which(nchar(annotations$annot) > 20)
    annotations$annot[who.huge] <- paste0(str_sub(temp$REF[who.huge],1,10),"...", " -> ",
                                          str_sub(temp$ALT[who.huge],1,10),"...")
    
    who.huge <- which(nchar(annotations$HGVS.c) > 25)
    annotations$HGVS.c[who.huge] <- paste0(str_sub(temp$HGVS.c[who.huge],1,25),"...")
    
    
    who.huge <- which(nchar(annotations$HGVS.p) > 25)
    annotations$HGVS.p[who.huge] <- paste0(str_sub(temp$HGVS.p[who.huge],1,25),"...")
    
  }
  
  annotations <- right_join(annotations, deletions[, c("match", "mutation")])
  
  # row.names(annotations) <- annotations$match
  # 
  # annotations <- annotations[, -c("match")]
  
  colnames(data_specific) <- timepoints$heat_labels
  
  return(list(annotations = annotations,
              data_specific = data_specific))
  
}



check_amino <- function(input_data,
                        others,
                        timepoints){
  
  input_data$POS <- as.numeric(input_data$POS)
  others$HGVS.p <- paste0("p.", others$ref_aa, others$codon_num, others$alt_aa)
  
  specific_data <- input_data[which((input_data$HGVS.p %in% others$HGVS.p) & 
                                    (input_data$Gene_Name %in% others$gene)), ]
  
  specific_data <- specific_data[order(specific_data$POS, decreasing = F), ]
  
  data_specific <- matrix(0, nrow = nrow(others), ncol = nrow(timepoints))
  row.names(data_specific) <- paste0(others$gene, "_", others$HGVS.p)
  colnames(data_specific) <- timepoints$samples
  
  #initialize annotations
  annotations <- data.table(annot = character(nrow(others)),
                            genes = character(nrow(others)),
                            HGVS.c = character(nrow(others)),
                            HGVS.p = character(nrow(others)),
                            POS = character(nrow(others)))
  
  row.names(annotations) <- paste0(others$gene, "_", others$HGVS.p)
  
  for (k in timepoints$samples){
    
    temp <- specific_data[which(specific_data$sample == k), 
                          c("POS", "REF", "ALT", "TYPE", "Gene_Name", "HGVS.c", "HGVS.p", "AF")]
    
    temp$length <- unlist(str_split(temp$HGVS.c, ">"))[2]
    temp$length <- nchar(temp$length)
    
    temp$dif <- stringdist(temp$REF, temp$ALT, method = "hamming")
    
    temp <- temp[which(temp$length == temp$dif), ]
    
    for  (i in c(1:nrow(others))){
      
      removed <- which(temp$HGVS.p == others$HGVS.p[i] & temp$Gene_Name != others$gene[i])
      if (length(removed) > 0)
      temp <- temp[-removed, ]
    
    }
    
    temp <- group_by(temp, HGVS.p, Gene_Name) %>%
      summarise(POS, REF,
                ALT = paste(ALT, collapse = ","), Gene_Name, HGVS.c, HGVS.p ,
                AF = sum(AF),) %>%
      as.data.table %>%
      unique()
    
    temp$label <- paste0(temp$Gene_Name,"_", temp$HGVS.p)
    data_specific[as.character(temp$label), k] <- temp$AF
    
    to.be.annotated <- which(row.names(annotations) %in% temp$label)
    temp$label <- reorder.factor(temp$label, new.order=row.names(annotations))
    new_temp <- temp %>%
       arrange(label)
    
    annotations[to.be.annotated, c("annot")] <- paste0(new_temp$REF, " -> ", new_temp$ALT)
    annotations$genes[to.be.annotated] <- new_temp$Gene_Name
    annotations$HGVS.c[to.be.annotated] <- new_temp$HGVS.c
    annotations$HGVS.p[to.be.annotated] <- new_temp$HGVS.p
    annotations$POS[to.be.annotated] <- new_temp$POS
    
    who.huge <- which(nchar(annotations$annot) > 20)
    annotations$annot[who.huge] <- paste0(str_sub(new_temp$REF[who.huge],1,10),"...", " -> ",
                                          str_sub(new_temp$ALT[who.huge],1,10),"...")
    
    who.huge <- which(nchar(annotations$HGVS.c) > 25)
    annotations$HGVS.c[who.huge] <- paste0(str_sub(new_temp$HGVS.c[who.huge],1,25),"...")
    
    
    who.huge <- which(nchar(annotations$HGVS.p) > 25)
    annotations$HGVS.p[who.huge] <- paste0(str_sub(new_temp$HGVS.p[who.huge],1,25),"...")
    
  }
  
  annotations <- right_join(annotations, others[, c("HGVS.p", "mutation")])
  
  annotations <- annotations[order(annotations$mutation), ]
  others <- others[order(others$mutation), ]
  row.names(annotations) <- paste0(others$gene, "_", others$HGVS.p)
  
  colnames(data_specific) <- timepoints$heat_labels
  
  return(list(annotations = annotations,
              data_specific = data_specific))
}



create_dotplot <- function(heat_matrix,
                           new_annotations,
                           use_nas,
                           folder,
                           specific_fp){
  
  if (use_nas == TRUE) {
    
    for(i in c(1:ncol(heat_matrix))){
      heat_matrix[(heat_matrix[,i] == 0), i] <- NA
    }
    
  }
  
  data_dot <- data.table(AF = as.vector(heat_matrix),
                         time = c(rep(0, nrow(heat_matrix)),
                                  rep(1, nrow(heat_matrix)),
                                  rep(2, nrow(heat_matrix)),
                                  rep(3, nrow(heat_matrix)),
                                  rep(4, nrow(heat_matrix)),
                                  rep(5, nrow(heat_matrix)),
                                  rep(6, nrow(heat_matrix))),
                         variant = rep(new_annotations$mutation, ncol(heat_matrix)))
  
 
  
  plot <- ggplot(data_dot, aes(x=time, y= AF, color = variant, group = variant)) + 
    geom_point(shape=19, size = 1, position = position_dodge(width = 0.1))+  
    geom_line(linetype = "dashed", alpha = 0.2)+
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x  =element_text( size = 11, color = "black"),
          axis.text.y = element_text( size = 16, color = "black"),
          axis.title.x = element_text( size = 18, color = "black"),
          axis.title.y = element_text( size = 18, color = "black")
          #legend.position = "none"
    ) +
    scale_x_continuous(labels=colnames(heat_matrix), breaks = c(0,1,2,3,4,5,6)) +
    geom_smooth(se = F,  formula = y ~ poly(x,4), method = lm, 
                color = "black", aes(x = time, y = AF), inherit.aes = FALSE)  
 
  save_image(print(plot),paste0(folder, "/dotplot_", specific_fp,".png"), 
             width = 14, height = 7,res = 200)
 
}