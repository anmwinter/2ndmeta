library(ggplot2)
library(ggdendro)
library(gridExtra)

setwd("~/Desktop/")
gc.size <- read.delim("~/Desktop/NOV 2013 GC vs Size Bacteria.csv")
head(gc.size)

p <- ggplot(gc.size,aes(x=Size..Mb., y=GC., color=Phyla)) + geom_point(shape=1) +
  ggtitle("GC Content versus Size of Genome for Novemeber 2013 Complete Bacterial Genomes") +
  xlab("Genome Size (Mb)")+
  ylab("GC Content")+
  guides(col = guide_legend(nrow = 12))+
  theme_bw() +
  theme(legend.position = "bottom")+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
  )

print(p)