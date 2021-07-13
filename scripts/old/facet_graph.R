library(ggplot2)
library(dplyr)

data.hist <- data.table(abs.total.dif = c(),
                        n = c())
maxT <- max(stats.table$max.abs.total.dif)

for (l in c(1:3)){
  temp <- data.table(abs.total.dif = seq(0,maxT,1),n = 0)
  
  one.table <- data.graphs[[l]]
  data <- count(one.table, abs.total.dif)
  
  who <- setdiff(c(seq(0,maxT,1)),data$abs.total.dif)
  
  
  data.hist <- rbind(data.hist,data)
  print(nrow(data))
  if (length(who) != 0){
    data.hist <- rbind(data.hist, data.table(abs.total.dif = who,
                                             n = 0))
  }
}

data.hist$comparison <- c(replicate(maxT+1,"s100-s200"),
                          replicate(maxT+1,"s100-s500"),
                          replicate(maxT+1,"s200-s500"))

p <- data.hist %>%
  ggplot( aes(x=abs.total.dif, y = (n/8794), fill= as.factor(abs.total.dif))) + ylim(c(0,0.4))+
  geom_histogram( stat = "identity", width = 0.8) +
  #scale_fill_brewer(palette = "Accent") +
  #ggtitle(paste0("Abs_Total_Dif_",names(data.graphs)[l])) +
  theme_bw() + 
  geom_text(label = data.hist$n,nudge_y = 0.015  ,size = 2.6, angle = 90)+
  theme(
    #plot.title = element_text(size=10),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, vjust = -0.05,hjust = 0.5),
    axis.title.y = element_text(size = 14, vjust = 1.5,hjust = 0.5)
  ) +
  xlab("Absolute difference values") + ylab("Frequency")+ 
  theme(legend.position="none") + scale_x_continuous(breaks = seq(0,maxT,3))+
  facet_grid(. ~ comparison, labeller=label_both)+
  theme(strip.text.x = element_text(size=13, color="black"),
        strip.background = element_rect(colour="black", fill="white", 
                                        size=2, linetype="solid")) # +
#scale_y_continuous(breaks = seq(0,1,0.2))
# scale_x_continuous(breaks = seq(1, max(data.hist$abs.total.dif),1))

saveImageHigh::save_image(print(p), file.name =paste0(outputs,"/facet_all.png"), width = 10, height = 6)
