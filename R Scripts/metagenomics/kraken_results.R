library(ggplot2)
library(plyr)
library(reshape2)

setwd("~/Documents/basespace/082815DN1genomes-26141121/kraken_data_report/")

filenames <- list.files()

read_csv_filename <- function(filename){
  ret <- read.csv(filename, header=FALSE, sep="\t")
  ret$sample_id <- filename
  ret
}

import.list <- llply(filenames, read_csv_filename)

combined <- do.call("rbind", import.list)

head(combined)

domain <- subset(combined, V4=="D")
phylum <- subset(combined, V4=="P")
class <- subset(combined, V4=="C")

p1 = "proteobacteria"
proteos <- subset(combined, grepl(p1, V6))
proteos

v1 = "viridae"
virus_order <- subset(combined, grepl(v1, V6))
virus_order

all_taxa <- rbind(phylum,proteos,virus_order)
all_taxa

row_sub = apply(all_taxa, 1, function(row) all(row !="Proteobacteria" ))
all_taxa <- all_taxa[row_sub,]

head(all_taxa)

taxa_cols <- c("sample_id","V6","V2")
taxa_for_plot <- all_taxa[taxa_cols]
taxa_for_plot

ggplot(taxa_for_plot)+ geom_bar(aes(x=V6,y=log(V2),fill=sample_id), stat="identity") +
  coord_flip() +
  theme_bw() +
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    plot.background = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(linetype=3, color="darkgray"),
    axis.text.y=element_text(size=rel(0.5)))

