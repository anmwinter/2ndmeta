# Take the hmmer hits file against the AntiResFam hmmer profile (http://www.dantaslab.org/resfams/) 
# and convert it to best hits, counts, and a human readable annotation.

# Reads in the file produced by:
# hmmsearch --cpu 6 --domtblout SRR681084_snottie_resfams.txt 
# ~/Documents/hmmer_profiles/antiresfams/Resfams-full.hmm SRR681084_snottite.faa  > /dev/null

hmmer_hits <- read.table("~/Desktop/foam_test/SRR681084_snottie_resfams.txt",sep="",
                   header=FALSE, skip=0)

library(plyr)
# Rename the column titles 
hmmer_hits <- rename(hmmer_hits, c(V1="sequence", V4="query name",V7="E-value",V8="score"))
head(hmmer_hits)

# Best practices for best hits put the hmmer score somewhere bewteen 20 and 25. Typically you see
# 20 or 25 in the literature. 
hmmer_best_hits <- subset(hmmer_hits, score>25)

# Create a temp file to hold our renamed best hits
hmmer_count.tmp <- hmmer_best_hits[, c("sequence", "query name")]

# Here we summarize the data to get our best hits counts
hmmer_BH.counts <- as.data.frame(summary(hmmer_count.tmp$`query name`))

# The data comes through as named integars so we wrangled it into a dataframe above.
# So we need to rename the columns. 
hmmer_BH.counts$query_id <- rownames(hmmer_BH.counts) 
hmmer_BH.counts <- rename(hmmer_BH.counts, c("summary(hmmer_count.tmp$`query name`)"="counts"))
hmmer_BH.counts



