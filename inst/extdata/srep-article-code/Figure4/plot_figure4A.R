


# citation("data.table")
# citation("stringr")
# citation("ggplot2")
# citation("ggsci")


library(data.table)
library(stringr)

rm(list = ls())

clinical_data = readxl::read_excel("Corrected_tables_v2.xlsx", sheet = "clinical")

clinical_data = as.data.table(clinical_data)

for(i in which(is.na(clinical_data$Period))) {
  
  clinical_data[i, ]$Period = clinical_data[i - 1, ]$Period
  clinical_data[i, ]$`No of patients / period` = clinical_data[i - 1, ]$`No of patients / period`
  
}

clinical_data$Period = factor(clinical_data$Period, levels = unique(clinical_data$Period))

clinical_data$Percentage = clinical_data$Percentage / 100

library(ggplot2)
library(ggsci)

clinical_data$label_position = clinical_data$Percentage + 0.06


gr1 = ggplot(data = clinical_data, aes(x = Percentage, y = Period, fill = `Lineages Detected`)) +
  
  geom_bar(position = position_dodge2(width = 1, preserve = "single"), 
           stat = "identity") +
    
    scale_fill_npg() +
  
  geom_text(aes(x = label_position, y = Period, label = `Lineages Detected`),
            position = position_dodge2(width = 1, preserve = "single"),
            size = 2.3,
            hjust = 0.5,
            vjust = 0.5) +
  
  scale_x_continuous(labels = scales::percent,
                     position = "bottom",
                     
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     
                     expand = c(0, 0, 0, 0.05)) +
  
  theme_minimal() +
  
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
        axis.ticks = element_line(),
        
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        
        legend.position = "bottom",
        
        axis.line = element_line())

