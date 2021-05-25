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
  deletions$match[who_one] <- paste0(as.numeric(deletions$codon_num[who_one]), "del")
  #deletions$match[who_one] <- paste0(as.numeric(deletions$codon_num[who_one]) + 1, "del")
  
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


pangolin_deletions <- function(input_data,
                               deletions,
                               timepoints,
                               aa_abbreviations){
  
  input_data$POS <- as.numeric(input_data$POS)
  deletions$position <- as.numeric(deletions$position)
  deletions$position <- deletions$position - 1
  
  # find matches
  specific_data <- input_data[which(input_data$POS %in% deletions$position), ]
  specific_data <- specific_data[which(specific_data$TYPE == "del"), ]
  
  specific_data <- specific_data[order(specific_data$POS, decreasing = F), ]
  
  data_specific <- matrix(0, nrow = nrow(deletions), ncol = nrow(timepoints))
  row.names(data_specific) <- deletions$position
  colnames(data_specific) <- timepoints$samples
  
  #initialize annotations
  annotations <- data.table(annot = character(nrow(deletions)),
                            genes = character(nrow(deletions)),
                            HGVS.c = character(nrow(deletions)),
                            HGVS.p = character(nrow(deletions)),
                            POS = character(nrow(deletions)))
  
  row.names(annotations) <- deletions$position
  
  for (k in timepoints$samples){
    
    temp <- specific_data[which(specific_data$sample == k), 
                          c("POS", "REF", "ALT", "TYPE", "Gene_Name", "HGVS.c", "HGVS.p", "AF")]
    
    #temp$length <- nchar(temp$REF) - nchar(temp$ALT) + 1
    
    who.one <- which(deletions$length == 1)
    one.replacement <- deletions$POS[who.one]
    not.really.ones <- which((temp$POS %in% one.replacement) & 
                                ((nchar(temp$REF))-(nchar(temp$ALT)) != 1))
     
    who.many <- which(deletions$length > 1)
    many.replacement <- deletions$position[who.many]
    not.really.manys <- which((temp$POS %in% many.replacement) & 
                                 ((nchar(temp$REF) - nchar(temp$ALT)) != deletions$length))
     
    removed <-c(not.really.ones, not.really.manys)
     
    if (length(removed > 0)){
      temp <- temp[-removed, ]
    }
    
    
    temp <- group_by(temp, POS) %>%
      summarise(POS, REF,
                ALT = paste(ALT, collapse = ","), Gene_Name, HGVS.c, HGVS.p ,
                AF = sum(AF)) %>%
      as.data.table %>%
      unique()
    
    data_specific[as.character(temp$POS), k] <- temp$AF
    
    to.be.annotated <- which(row.names(annotations) %in% as.character(temp$POS))
    
    data_specific[to.be.annotated, k] <- temp$AF
    
    annotations[to.be.annotated, c("annot")] <- paste0(temp$REF, " -> ", temp$ALT)
    annotations$genes[to.be.annotated] <- temp$Gene_Name
    annotations$HGVS.c[to.be.annotated] <- temp$HGVS.c
    annotations$HGVS.p[to.be.annotated] <- temp$HGVS.p
    annotations$POS[to.be.annotated] <- temp$POS
    
    who.huge <- which(nchar(annotations$annot) > 20)
    annotations$annot[who.huge] <- paste0(str_sub(temp$REF[who.huge],1,10),"...", " -> ",
                                          str_sub(temp$ALT[who.huge],1,10),"...")
    
    who.huge <- which(nchar(annotations$HGVS.c) > 25)
    annotations$HGVS.c[who.huge] <- paste0(str_sub(temp$HGVS.c[who.huge],1,25),"...")
    
    
    who.huge <- which(nchar(annotations$HGVS.p) > 25)
    annotations$HGVS.p[who.huge] <- paste0(str_sub(temp$HGVS.p[who.huge],1,25),"...")
    
  }
  
  deletions$position <- as.character(deletions$position)
  
  annotations <- right_join(annotations, deletions[, c("label", "position")],
                            by = c("POS" = "position"))
  
  deletions <- deletions[order(deletions$label), ]
  annotations$label <- reorder.factor(annotations$label, new.order = deletions$label)
  
  new_annotations <- annotations %>%
    arrange(label)
  annotations <- new_annotations
  
  row.names(annotations) <- deletions$position
  
  colnames(data_specific) <- timepoints$heat_labels
  
  annotations$new <- row.names(annotations)
  annotations$new <- reorder.factor(annotations$new, new.order=row.names(data_specific))
  new_annotations <- annotations %>%
    arrange(new)
  new_annotations <- new_annotations[, -c("new")]
  
  row.names(new_annotations) <- new_annotations$label
  
  return(list(annotations = new_annotations,
              data_specific = data_specific))
  
}


veo_deletions <- function(input_data,
                          deletions,
                          timepoints,
                          aa_abbreviations){
  
  deletions$HGVS.p <- paste0("p.", deletions$HGVS.p)
  
  # find matches
  specific_data <- input_data[which((input_data$HGVS.p %in% deletions$HGVS.p) &
                                (input_data$Gene_Name %in% deletions$gene)), ]
  
  data_specific <- matrix(0, nrow = nrow(deletions), ncol = nrow(timepoints))
  row.names(data_specific) <- deletions$HGVS.p
  colnames(data_specific) <- timepoints$samples
  
  #initialize annotations
  annotations <- data.table(annot = character(nrow(deletions)),
                            genes = character(nrow(deletions)),
                            HGVS.c = character(nrow(deletions)),
                            HGVS.p = character(nrow(deletions)),
                            POS = character(nrow(deletions)))
  
  row.names(annotations) <- deletions$HGVS.p
  
  for (k in timepoints$samples){
    
    temp <- specific_data[which(specific_data$sample == k), 
                          c("POS", "REF", "ALT", "TYPE", "Gene_Name", "HGVS.c", "HGVS.p", "AF")]
    
    temp <- group_by(temp, POS) %>%
      summarise(POS, REF,
                ALT = paste(ALT, collapse = ","), Gene_Name, HGVS.c, HGVS.p ,
                AF = sum(AF)) %>%
      as.data.table %>%
      unique()
    
    data_specific[as.character(temp$HGVS.p), k] <- temp$AF
    
    to.be.annotated <- which(row.names(annotations) %in% as.character(temp$HGVS.p))
    
    data_specific[to.be.annotated, k] <- temp$AF
    
    annotations[to.be.annotated, c("annot")] <- paste0(temp$REF, " -> ", temp$ALT)
    annotations$genes[to.be.annotated] <- temp$Gene_Name
    annotations$HGVS.c[to.be.annotated] <- temp$HGVS.c
    annotations$HGVS.p[to.be.annotated] <- temp$HGVS.p
    annotations$POS[to.be.annotated] <- temp$POS
    
    who.huge <- which(nchar(annotations$annot) > 20)
    annotations$annot[who.huge] <- paste0(str_sub(temp$REF[who.huge],1,10),"...", " -> ",
                                          str_sub(temp$ALT[who.huge],1,10),"...")
    
    who.huge <- which(nchar(annotations$HGVS.c) > 25)
    annotations$HGVS.c[who.huge] <- paste0(str_sub(temp$HGVS.c[who.huge],1,25),"...")
    
    
    who.huge <- which(nchar(annotations$HGVS.p) > 25)
    annotations$HGVS.p[who.huge] <- paste0(str_sub(temp$HGVS.p[who.huge],1,25),"...")
    
  }
  
  #deletions$position <- as.character(deletions$position)
  
  annotations <- right_join(annotations, deletions[, c("label", "HGVS.p")],
                            by = c("HGVS.p" = "HGVS.p"))
  
  deletions <- deletions[order(deletions$label), ]
  annotations$label <- reorder.factor(annotations$label, new.order = deletions$label)
  
  new_annotations <- annotations %>%
    arrange(label)
  annotations <- new_annotations
  
  row.names(annotations) <- deletions$HGVS.p
  
  colnames(data_specific) <- timepoints$heat_labels
  
  annotations$new <- row.names(annotations)
  annotations$new <- reorder.factor(annotations$new, new.order=row.names(data_specific))
  new_annotations <- annotations %>%
    arrange(new)
  new_annotations <- new_annotations[, -c("new")]
  
  
  row.names(new_annotations) <- new_annotations$label
  
  return(list(annotations = new_annotations,
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
    
    
    who_more <- which(str_detect(temp$HGVS.c, "ins"))
    
    temp$length[who_more] <- unlist(str_split(temp$HGVS.c[who_more], "ins"))[2]
    
    temp$length[who_more] <- nchar(temp$length[who_more])
    
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


pangolin_amino <- function(input_data,
                        others,
                        timepoints){
  
  input_data$POS <- as.numeric(input_data$POS)
  others$HGVS.p <- paste0("p.", others$HGVS.p)
  
  input_data$Gene_Name <- toupper(input_data$Gene_Name)
  others$gene <- toupper(others$gene)
  input_data$match <- paste0(input_data$Gene_Name, "_", input_data$HGVS.p)
  others$match <- paste0(others$gene, "_", others$HGVS.p)
  
  specific_data <- input_data[which(input_data$match %in% others$match), ] 
  
  specific_data <- specific_data[, -length(specific_data)]
  input_data <- input_data[, -length(input_data)]
  others <- others[, -length(others)]
  
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
    
    who_more <- which(str_detect(temp$HGVS.c, "ins"))
    
    temp$length[who_more] <- unlist(str_split(temp$HGVS.c[who_more], "ins"))[2]
    
    temp$length[who_more] <- nchar(temp$length[who_more])
    
    temp$dif <- stringdist(temp$REF, temp$ALT, method = "hamming")
     
    temp <- temp[which(temp$length == temp$dif), ]
    
    # 
    # for  (i in c(1:nrow(others))){
    #   
    #   removed <- which(temp$HGVS.p == others$HGVS.p[i] & temp$Gene_Name != others$gene[i])
    #   if (length(removed) > 0)
    #     temp <- temp[-removed, ]
    #   
    # }
    
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
  
  annotations <- right_join(annotations, others[, c("HGVS.p", "label")])
  
  # annotations <- annotations[order(annotations$label), ]
  others <- others[order(others$label), ]
  
  annotations$label <- reorder.factor(annotations$label, new.order = others$label)
  new_annotations <- annotations %>%
    arrange(label)
  annotations <- new_annotations
  
  row.names(annotations) <- paste0(others$gene, "_", others$HGVS.p)
  
  #colnames(data_specific) <- timepoints$heat_labels
  
  #row.names(annotations) <- annotations$label
  
  colnames(data_specific) <- timepoints$heat_labels
  
  annotations$new <- row.names(annotations)
  annotations$new <- reorder.factor(annotations$new, new.order=row.names(data_specific))
  new_annotations <- annotations %>%
    arrange(new)
  new_annotations <- new_annotations[, -c("new")]
  
  row.names(new_annotations) <- new_annotations$label
  
  return(list(annotations = new_annotations,
              data_specific = data_specific))
}



create_dotplot <- function(heat_matrix,
                           new_annotations,
                           use_nas,
                           plots_folder,
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
                                  rep(6, nrow(heat_matrix)),
                                  rep(7, nrow(heat_matrix)),
                                  rep(8, nrow(heat_matrix)),
                                  rep(9, nrow(heat_matrix)),
                                  rep(10, nrow(heat_matrix)),
                                  rep(11, nrow(heat_matrix)),
                                  rep(12, nrow(heat_matrix)),
                                  rep(13, nrow(heat_matrix))),
                         variant = rep(new_annotations$mutation, ncol(heat_matrix)))
  
 
  line_labels <- data.table(label = round(colMeans(heat_matrix, na.rm = T), 4),
                            x = seq(0, ncol(heat_matrix)-1, 1),
                            y = colMeans(heat_matrix, na.rm = T))
  
  line_labels <- rbind(line_labels, data.table(label = rep(" ", nrow(data_dot) - ncol(heat_matrix)),
                                               x = NA,
                                               y = NA))
  
  plot <- ggplot(data = data_dot, aes(x=time, y= AF, color = data_dot$variant, group = data_dot$variant)) + 
    geom_point(shape=19, size = 1, position = position_dodge(width = 0.1))+  
    geom_line(linetype = "dashed", alpha = 0.2)+
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x  =element_text( size = 13, color = "black", angle = 90),
          axis.text.y = element_text( size = 16, color = "black"),
          axis.title.x = element_text( size = 18, color = "black"),
          axis.title.y = element_text( size = 18, color = "black")
          #legend.position = "none"
    ) +
    guides(fill = guide_legend(title = "variant")) +
    scale_x_continuous(labels=colnames(heat_matrix), breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)) +
    geom_smooth(se = F,  formula = y ~ poly(x,4), method = lm, 
                color = "black", aes(x = time, y = AF), inherit.aes = FALSE) +
    geom_label_repel(data = line_labels, aes (x = x, y = y, label = label), 
                     color = "black", size = 3) +
    labs(color = "variant")
    
  dir.create(paste0(plots_folder, "/dotplots"))
  
  if (use_nas == TRUE){
    save_image(print(plot), paste0(plots_folder, "/dotplots/dotplot_nas_", specific_fp,".png"), 
               width = 18, height = 7,res = 200)
  } else {
    save_image(print(plot), paste0(plots_folder, "/dotplots/dotplot_", specific_fp,".png"), 
               width = 18, height = 7,res = 200)
  }
  
 
}


create_boxplot <- function(heat_matrix,
                           new_annotations,
                           plots_folder,
                           specific_fp){
  
  data_box <- data.table(AF = as.vector(heat_matrix),
                         time = c(rep("0", nrow(heat_matrix)),
                                  rep("1", nrow(heat_matrix)),
                                  rep("2", nrow(heat_matrix)),
                                  rep("3", nrow(heat_matrix)),
                                  rep("4", nrow(heat_matrix)),
                                  rep("5", nrow(heat_matrix)),
                                  rep("6", nrow(heat_matrix)),
                                  rep("7", nrow(heat_matrix)),
                                  rep("8", nrow(heat_matrix)),
                                  rep("9", nrow(heat_matrix)),
                                  rep("91", nrow(heat_matrix)),
                                  rep("92", nrow(heat_matrix)),
                                  rep("93", nrow(heat_matrix)),
                                  rep("94", nrow(heat_matrix))))
                         # time = c(rep(colnames(heat_matrix)[1], nrow(heat_matrix)),
                         #          rep(colnames(heat_matrix)[2], nrow(heat_matrix)),
                         #          rep(colnames(heat_matrix)[3], nrow(heat_matrix)),
                         #          rep(colnames(heat_matrix)[4], nrow(heat_matrix)),
                         #          rep(colnames(heat_matrix)[5], nrow(heat_matrix)),
                         #          rep(colnames(heat_matrix)[6], nrow(heat_matrix)),
                         #          rep(colnames(heat_matrix)[7], nrow(heat_matrix))))
  
  labels <- data.table(names = paste0(colSums(heat_matrix != 0), "/", nrow(heat_matrix)),
                       pos_y = apply(heat_matrix, 2, quantile)[4, ],
                       pos_x = seq(0,13,1))
  
  plot <- ggplot(data = data_box , aes(x=time, y=AF, fill=time)) +
    geom_boxplot(position=position_dodge(width=0.8), width = 0.4) +
    #geom_jitter(color="black", size=1, alpha=1) +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme_bw() + 
    theme(
      legend.position="none",
      plot.title = element_text(size=18),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x  =element_text( size = 13, color = "black", angle = 90),
      axis.text.y = element_text( size = 16, color = "black"),
      axis.title.x = element_text( size = 18, color = "black"),
      axis.title.y = element_text( size = 18, color = "black")
    ) +
    ggtitle((toupper(specific_fp)))  +
    xlab("Time") +
    geom_text(data=labels, aes(label=names, y = pos_y + 0.1, x = pos_x + 1.25), 
              size=4.5, position=position_dodge(width=0.1), inherit.aes = FALSE) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2), 
                       labels = c("0", "0.2", "0.4", "0.6", "0.8", "1", " ")) +
    scale_x_discrete(labels=colnames(heat_matrix), 
                     breaks = c("0","1","2","3","4","5","6", "7", "8", "9", "91", "92", "93", "94")) 

  
  dir.create(paste0(plots_folder, "/boxplots"))
  save_image(print(plot), paste0(plots_folder, "/boxplots/boxplot_", specific_fp,".png"), 
             width = 14, height = 7,res = 200)
  
  
}