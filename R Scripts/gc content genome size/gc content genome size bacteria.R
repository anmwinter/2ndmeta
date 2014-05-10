library(ggplot2)
library(gridExtra)
library(vegan)
library(FactoMineR)
library(ggdendro)
library(pgirmess)
library(BayesFactor)
library(BayesianFirstAid)

setwd("~/Desktop/PhD_Project/genome level bacteria stats/April 2014/")
gc_genome <- read.delim("~/Desktop/PhD_Project/genome level bacteria stats/April 2014/prokaryotes April 2014.csv",na.strings="-" ,sep="\t")
head(gc_genome)
tail(gc_genome)
summary(gc_genome)

filtered_gc_genome <- gc_genome[,c(1,5,6,7,8,19)]
head(filtered_gc_genome)
summary(filtered_gc_genome)
filtered_gc_genome <- na.omit(filtered_gc_genome)
summary(filtered_gc_genome)

complete_genomes <- filtered_gc_genome[filtered_gc_genome$Status=="Complete",]
summary(complete_genomes)

write.table(complete_genomes, file = "summary genome stats.csv", sep = ";",quote=FALSE)

#Global GC Content, Genome Size Bayesian Correlation

qplot(GC., Size..Mb., data = complete_genomes, shape = I(1)) #+ facet_grid(Geotype ~ ., scales = "free")

size_vs_gc <- ggplot(complete_genomes,aes(x=Size..Mb., y=GC.,)) + geom_point(shape=1) +
  ggtitle("GC Content vs Genzome Size - Bacterial and Euryarchaeota Genomes from NCBI") +
  xlab("Genzome Size (Mb)")+
  ylab("% GC Content")+
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
  )
size_vs_gc

#Bayesian First Aid Pearson's Correlation Coefficient Test
#on completed genomes

fit_meta <- bayes.cor.test( ~ Size..Mb. + GC., data = complete_genomes)
fit_meta
plot(fit_meta)
#Shows you the distribution we would expect a new data point to have
#Light blue oval is the 95% density region
#Dark blue oval is the 50% density region
