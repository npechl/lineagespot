

rm(list = ls())


# citation("data.table")
# citation("stringr")
# citation("xlsx")
# citation("ggplot2")
# citation("dplyr")
# citation("hrbrthemes")
# citation("ggsci")
# citation("extrafont")
# citation("patchwork")



library(data.table)
library(stringr)
library(xlsx)


# read inputs --------------
meta_file <- "meta-data.xlsx"
meta_data <- xlsx::read.xlsx(meta_file, sheetIndex = 1)[1:14, ]

timepoints <- data.table(samples = meta_data$Sample_ID,
                         heat_labels = meta_data$Time.length_.sampling.)


rm(meta_data, meta_file)

unique_mutations <- read.csv("unique_mutations_B.1.1.7.csv")
all_mutations <- fread("outbreakinfo_mutation_report_data_B.1.1.7.tsv")[, 1]
hits <- read.csv("results_B.1.1.7.csv")

# create heatmap --------------
hits$name <- paste(hits$Gene_Name, hits$HGVS.p, sep = ":")
all_new_mu <- unique(hits$name)

entire_table <- matrix(0,
                       nrow = length(all_new_mu),
                       ncol = nrow(timepoints))


colnames(entire_table) <- timepoints$samples
row.names(entire_table) <- all_new_mu

for (k in timepoints$samples) {
  
  temp <- hits[which(hits$sample == k), ]
  row.names(temp) <- temp$name
  entire_table[temp$name, k] <- temp$AF
  
}

colnames(entire_table) <- timepoints$heat_labels

#entire_table$mutation <- row.names(entire_table)

# create dotplot --------
library(ggsci)

colors = c(pal_igv("default")(51)) #pal_ucscgb("default")(26)
#pal_d3("category20")(20)

data_dot <- data.table(AF = as.vector(entire_table),
                       time = c(rep(0, nrow(entire_table)),
                                rep(1, nrow(entire_table)),
                                rep(2, nrow(entire_table)),
                                rep(3, nrow(entire_table)),
                                rep(4, nrow(entire_table)),
                                rep(5, nrow(entire_table)),
                                rep(6, nrow(entire_table)),
                                rep(7, nrow(entire_table)),
                                rep(8, nrow(entire_table)),
                                rep(9, nrow(entire_table)),
                                rep(10, nrow(entire_table)),
                                rep(11, nrow(entire_table)),
                                rep(12, nrow(entire_table)),
                                rep(13, nrow(entire_table))),
                       variant = rep(row.names(entire_table), ncol(entire_table)))
# 
# 
# line_labels <- data.table(label = round(colMeans(heat_matrix, na.rm = T), 4),
#                           x = seq(0, ncol(heat_matrix)-1, 1),
#                           y = colMeans(heat_matrix, na.rm = T))
# 
# line_labels <- rbind(line_labels, data.table(label = rep(" ", nrow(data_dot) - ncol(heat_matrix)),
#                                              x = NA,
#                                              y = NA))

library(ggplot2)

my.colors = pal_material("grey", 
                         n = 2 * length(unique(data_dot$variant)),
                         reverse = TRUE)(2 * length(unique(data_dot$variant)))

my.colors = my.colors[seq(1, 2 * length(unique(data_dot$variant)), by = 2)]


plot = ggplot(data = data_dot, 
       aes(x=time, y= AF, 
           # color = data_dot$variant, 
           group = data_dot$variant)) + 
    
  geom_point(size = 3, position = position_dodge(width = 0.3), 
             alpha = 0.2)+ 
    
    # scale_color_manual(values = my.colors) +
    
  ylab("Allele frequency") +
  # geom_line(linetype = "dashed", alpha = 0.2)+
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        
        axis.text.x  =element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y = element_text( size = 10),
        
        axis.title.x = element_blank(),
        axis.ticks = element_line(),
        axis.title.y = element_text( size = 11),
        # legend.title = element_text(size = 13),
        # legend.title.align = "top",
        legend.text = element_text(size = 10),
        
        axis.line = element_line(),
        
        legend.position = "bottom"
  ) +
    
  # guides(fill = guide_legend(title = "Mutation")) +
    
  scale_x_continuous(labels = colnames(entire_table), 
                     breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)) +
    
    scale_y_continuous(expand = c(0, 0, 0, 0.05)) +
    
  geom_smooth(se = TRUE,  
              formula = y ~ poly(x,4), 
              method = lm, 
              color = "#3C5488FF", 
              size = 2,
              aes(x = time, y = AF), 
              inherit.aes = FALSE)
    
  # scale_fill_manual(values=colors) +
  #  geom_label_repel(data = line_labels, aes (x = x, y = y, label = label), 
  #                   color = "black", size = 3) +
  # labs(color = "variant")

plot

# save_image(print(plot), paste0("tests_graphs/dotplot-3.png"), 
#            width = 17, height = 7,res = 200) 
# 
# save_as_pdf(print(plot),
#             paste0("tests_graphs/dotplot-3.pdf"), 
#             width = 17, height = 7)

# Try reversed dotplot -------------
# 
# data_dot2 <- data.table(AF = as.vector(entire_table),
#                        time = rep(timepoints$heat_labels, each = nrow(entire_table)),
#                        variant = rep(row.names(entire_table), ncol(entire_table)))
# # 
# 
# 
# plot2 <- ggplot(data = data_dot2, aes(x=data_dot$variant, y= AF, fill = time, group = time)) + 
#   geom_point(shape=21, size = 5, position = position_dodge(width = 0.3), alpha = 0.6)+ 
#   ylab("Allele Frequency") +
#   # geom_line(linetype = "dashed", alpha = 0.2)+
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x  =element_text( size = 13, color = "black", angle = 45, hjust = 1, vjust = 1),
#         axis.text.y = element_text( size = 16, color = "black"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text( size = 18, color = "black"),
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 13),
#         axis.line = element_line(colour = "black", 
#                                  size = 1,
#                                  linetype = "solid"),
#         legend.position = "bottom"
#   ) +
#   guides(fill = guide_legend(title = "Timeperiod")) +
#  # scale_x_continuous(labels=colnames(entire_table), breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)) +
#  # geom_smooth(se = F,  formula = y ~ poly(x,4), method = lm, 
# #              color = "black", aes(x = time, y = AF), inherit.aes = FALSE) +
#   scale_fill_manual(values=colors) +
#   #  geom_label_repel(data = line_labels, aes (x = x, y = y, label = label), 
#   #                   color = "black", size = 3) +
#   labs(color = "time")
# 
# plot2
# 
# save_image(print(plot2), paste0("tests_graphs/dotplot-4.png"), 
#            width = 17, height = 7,res = 200)  
# 
# save_as_pdf(print(plot2),
#             paste0("tests_graphs/dotplot-4.pdf"), 
#             width = 17, height = 7)
# Calculate vectors ----------

# Read amino acid abbreviations --------------

amino_acids <- "Amino_acid_Abbreviations.txt"
aa_abbreviations <- read.table(amino_acids, header = T, sep = ",")
rm(amino_acids)


# mean all 

mean_all <- colMeans(entire_table)

for (i in c(1:nrow(aa_abbreviations))){
  
  row.names(entire_table) <- str_replace_all( row.names(entire_table), 
                       aa_abbreviations$Three_Letter[i], 
                       aa_abbreviations$One_Letter[i])
}

row.names(entire_table) <- tolower( row.names(entire_table))

who <- which(str_detect(unique_mutations$x, "orf1a:|orf1b:"))
unique_mutations$x <- str_replace_all(unique_mutations$x,
                                      "orf1a:|orf1b:",
                                      "orf1ab:")

who <- which(str_detect(row.names(entire_table), "del") & 
               str_detect(row.names(entire_table), "_"))

to_be_changed <- row.names(entire_table)[who]  #gsub("[^0-9.-]", "",row.names(entire_table)[who])

to_be_changed <- str_split(to_be_changed, ":|_", simplify = T) 

to_be_changed[, 2] <- gsub("[^0-9.-]", "",to_be_changed[, 2]) 
to_be_changed[, 3] <- gsub("[^0-9.-]", "",to_be_changed[, 3]) 

row.names(entire_table)[who] <- paste0(to_be_changed[, 1], ":del", 
                                       to_be_changed[, 2], "/",
                                       to_be_changed[, 3])

table_unique <- entire_table[unique_mutations$x, ]

mean_unique <- colMeans(table_unique)

table_unique[table_unique == 0] <- 2

min_unique <- apply(table_unique, 2, min)
min_unique[min_unique == 2] <- 0

# write.csv(data.frame(lapply(mean_all, type.convert), stringsAsFactors=FALSE),
#           "tests_graphs/mean_all_B.1.1.7.csv",
#           row.names = F)
# 
# write.csv(data.frame(lapply(mean_unique, type.convert), stringsAsFactors=FALSE),
#           "tests_graphs/mean_unique_B.1.1.7.csv",
#           row.names = F)
# 
# write.csv(data.frame(lapply(min_unique, type.convert), stringsAsFactors=FALSE),
#           "tests_graphs/min_unique_B.1.1.7.csv",
#           row.names = F)

# Create lineplot -------------

# Libraries
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(ggsci)
library(extrafont)

loadfonts()

theme_set(theme_bw(base_size = 12, base_family = "Helvetica") + 
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

data <- data.table::data.table(timeperiods = c("2 - 14 December 2020",
                                               "5 - 11 February 2021",
                                               "12 - 18 February 2021",
                                               "19 - 25 February 2021",
                                               "26 February - 4 March 2021",
                                               "5 - 11 March 2021",
                                               "12 - 18 March 2021",
                                               "19 - 25 March 2021",
                                               "26 March - 1 April 2021",
                                               "2 - 8 April 2021",
                                               "9 - 15 April 2021",
                                               "16 - 22 April 2021",
                                               "23 - 29 April 2021",
                                               "30 April - 6 May 2021"),
                               avg_all_mutations = mean_all *100,
                               avg_unique_mutations = mean_unique *100,
                               min_unique_mutations = min_unique * 100,
                               clinical_data = c(0.00, 43.9, 53.1,
                                                 81, 83.2, 100,
                                                 100, 100, 100,
                                                 100, 100, 96.9,
                                                 98.5, 100))

data$timeperiods <- factor(data$timeperiods, levels =  unique(data$timeperiods))

data3 <- data.table::data.table(timepoints = rep(data$timeperiods, 2),
                                percentages = c(data$min_unique_mutations, data$clinical_data),
                                Type = c(rep("Min AF of existing unique mutations of wastewater data", nrow(data)),
                                         rep("Clinical data", nrow(data))))
data3$timepoints <- factor(data3$timepoints, levels =  unique(data3$timepoints))

data3$percentages = data3$percentages / 100

plot3 <- ggplot(data = data3 ,
                aes(x = timepoints, y = percentages, group = Type, color = Type)) +
  geom_line(size = 2, alpha = 1, aes(linetype = Type)) +
  geom_point(size = 4, position = position_dodge(width = 0.1),
             alpha = 1) +
  scale_color_manual(values = c("gray50", "gray70")) +
  ggtitle("") +
  ylab("Allele frequency") +
  xlab(" ") +
  theme_minimal() +
  # theme_set(theme_bw(base_family = "Arial Narrow")) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x  =element_text(angle = 45, vjust=1, hjust=1, size = 10),
        axis.text.y = element_text( size = 10),
        
        # axis.title.x = element_text( size = 18, color = "black", hjust =0.5),
        axis.title.y = element_text(hjust = 0.5, vjust = 1, size = 11),
        axis.line = element_line(),
        
        axis.ticks = element_line(),
        
        # legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 10),
        legend.position = "top") +
  
  # scale_y_continuous(breaks = seq(0, 100, 20), labels = paste0(seq(0, 100 /100, 20 /100))) +
  scale_x_discrete(breaks = unique(data3$timepoints), 
                   labels = unique(data3$timepoints)) 


library(patchwork)

plot3 / plot

