library(readr)
library(data.table)
library(stringr)

#the path of the folder with all the files to be compared
input.folder <- "Data/L1"

files <- list.files(path = input.folder,pattern = "tsv")

#create output folder
outputs <- "Outputs_L1"
dir.create(outputs)

files.list <- list()

total <- read_tsv(paste0(input.folder,"/",files[1]))[,1:3]
total  <-  total[order(total$Lineage, decreasing = F),]

#read and combine files
for (i in c(1:length(files))){
  
  one.file <- read_tsv(paste0(input.folder,"/",files[i]))
  files.list[[i]] <- one.file
  
  one.file <- one.file[,c("Tree.Overlap", "Total.Overlap" )]
  colnames(one.file) <- paste0(colnames(one.file),"_",files[i])
  
  total <- cbind(total, one.file)
  
}

num.files <- length(files)
colnames(total) <- str_remove_all(colnames(total),".tsv")

iter <- seq(4,(2*num.files+2),2)

#compare files
stats.table <- data.table(files = character(),
                          num.differences = numeric(),
                          max.abs.tree.dif = numeric(),
                          max.abs.total.dif = numeric())
#data.graphs <- list()
#names.list <- c()
#k <- 1
for (i in iter[1:length(iter)-1]){
  for (j in seq((i+2),(2*num.files+2),2)){
   
    one.table <- total[,c(1:3,i, i+1, j, j+1)]
    
    filename1 <- str_remove(colnames(total)[i],"Tree.Overlap_")
    filename2 <- str_remove(colnames(total)[j],"Tree.Overlap_")
    tree.dif.name <- paste0("Dif.Tree.Overlap_",filename1,"_",filename2)
    total.dif.name <- paste0("Dif.Total.Overlap_",filename1,"_",filename2)
    abs.tree.name <- paste0("Abs.",tree.dif.name)
    abs.total.name <- paste0("Abs.",total.dif.name)
    
    #tree.dif <- total[,i]-total[,j]
    abs.tree.dif <- abs(total[,i]-total[,j])
    #total.dif <- total[,i+1]-total[,j+1]
    abs.total.dif <- abs(total[,i+1]-total[,j+1])
    
    one.table <- cbind(one.table,abs.tree.dif ,abs.total.dif )
    #data.graphs[[k]] <- one.table
    #names.list <- c(names.list,paste0(filename1,"_",filename2))
    #k <-k+1
    one.table <- one.table[which(one.table[,9] != 0),]
    colnames(one.table) <- c(colnames(total)[c(1:3,i, i+1, j, j+1)],
                             abs.tree.name ,abs.total.name)
    
    one.table <- one.table[order(one.table[,9],decreasing = T),]
    
    if (nrow(one.table)>0) {
      
      write.table(one.table, paste0(outputs,"/differences_",filename1,"_",filename2,".csv"), sep = "\t", row.names = F)
      
    }
    
    stats.table <- rbind(stats.table, data.table(files = paste0(filename1,"_", filename2),
                                                 num.differences = nrow(one.table),
                                                 max.abs.tree.dif = max(one.table[,8]),
                                                 max.abs.total.dif = max(one.table[,9])))
    
  }
  
}

write.table(stats.table,paste0(outputs,"/info for all the diffeneces.csv"), sep = "\t", row.names = F)
#names(data.graphs) <- names.list

#names(files.list) <- files

#table <- matrix(nrow = length(files), ncol = length(files))
#dimnames(table) <- list(files, files)

#for (i in c(1:ncol(table))){
#  for (j in c(i:ncol(table))){
#    
#    table[j,i] <- identical(files.list[[i]][,c("Tree.Overlap", "Total.Overlap")], files.list[[j]][,c("Tree.Overlap", "Total.Overlap")])
#    
#  }
#}

#write.table(table, paste0(outputs,"/check_identical_files.txt"), sep = "\t", row.names = T)

