# Take the hmmer hits file against the EGGNOG databases and convert it to best hits, counts, and
# human readable annotation

hmmer_hits <- read.table("~/Desktop/foam_test/SRR681084_snottie_retorvir.txt",sep="",
                   header=FALSE, skip=0)

library(plyr)
hmmer_hits <- rename(hmmer_hits, c(V1="sequence", V4="query name",V7="E-value",V8="score"))
head(hmmer_hits)

hmmer_best_hits <- subset(hmmer_hits, score>25)

hmmer_count.tmp <- hmmer_best_hits[, c("sequence", "query name")]

hmmer_BH.counts <- as.data.frame(summary(hmmer_count.tmp$`query name`))

hmmer_BH.counts$query_id <- rownames(hmmer_BH.counts) 
hmmer_BH.counts <- rename(hmmer_BH.counts, c("summary(hmmer_count.tmp$`query name`)"="counts"))
hmmer_BH.counts

# Convert the eggnog query id to human readable

# Reads in the eggnog annotation file (a .tsv) and links the query id to the human readable annotation
eggnog_ids <- read.csv("~/Desktop/foam_test/Retroviridae.annotations.tsv",sep="\t",
                         header=TRUE, skip=0)

eggnog_ids <- rename(eggnog_ids, c("OG"="query_id"))

besthits_annotation <- merge(hmmer_BH.counts,eggnog_ids,by="query_id")

besthits_annotation


