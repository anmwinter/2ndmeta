library(tm)
library(SnowballC)
library(ggplot2)
library(ggdendro)
library(kernlab)
require(FactoMineR)
library(slam)
library(corrplot)
library(pvclust)
library(ape)
library(cluster)
library(topicmodels)
library(lsa)


# set working directory, ie. location of JSTOR DfR CSV
# files on the computer
setwd("~/Desktop/ICLSJSTOR_design_by_disipline/design_art/wordcounts")

# create a list of all the CSV files
myFiles <- list.files(pattern="*.csv|CSV")

# read in all the CSV files to an R data object
myData <-  lapply(myFiles, read.csv)

# assign file names to each dataframe in the list
names(myData) <- myFiles

# Here's the step where we turn the JSTOR DfR 'wordcount' into
# the 'bag of words' that's typically needed for topic modelling
# The R process is 'untable-ing' each CSV file into a
# list of data frames, one data frame per file
myUntabledData <- sapply(1:length(myData),
                         function(x) {rep(myData[[x]]$WORDCOUNTS, times = myData[[x]]$WEIGHT)})

names(myUntabledData) <- myFiles
sapply(myFiles,
       function (x) write.table(myUntabledData[x], file=paste(x, "txt", sep="."),
                                quote = FALSE, row.names = FALSE, eol = " " ))



rpg_corpus  <-Corpus(DirSource("~/Desktop/ICLSJSTOR_design_by_disipline/design_art/texts"),  readerControl = list(reader = readPlain, load=TRUE))
print(rpg_corpus)
summary(rpg_corpus)

rpg_corpus <- tm_map(rpg_corpus, removeNumbers) 
rpg_corpus <- tm_map(rpg_corpus, removePunctuation)
rpg_corpus <- tm_map(rpg_corpus, stripWhitespace) # Removes whitespace in the texts
rpg_corpus <- tm_map(rpg_corpus, tolower) # Makes all words lowercase
rpg_corpus <- tm_map(rpg_corpus, stemDocument) # Stems the document to reduce size
rpg_corpus <- tm_map(rpg_corpus, removeWords, stopwords("english"))
rpg_corpus <- tm_map(rpg_corpus, removeWords, stopwords("SMART"))
###What follows are the high frequecny words that don't matter
rpg_corpus <- tm_map(rpg_corpus, removeWords, ("time"))
rpg_corpus <- tm_map(rpg_corpus, removeWords, ("result"))
rpg_corpus <- tm_map(rpg_corpus, removeWords, ("year"))
rpg_corpus <- tm_map(rpg_corpus, removeWords, ("number"))
rpg_corpus <- tm_map(rpg_corpus, removeWords, ("american"))
rpg_corpus <- tm_map(rpg_corpus, removeWords, ("john"))


tdm <- DocumentTermMatrix(rpg_corpus) #Takes the corpus and makes it into a dtm
dim(tdm) #This shows the dimensions of the matrix. The second number is the number of words.
tdm <- removeSparseTerms(tdm, .50) #Removes words that occur in only a few documents. 
dim(tdm) #This shows the dimensions of the matrix. The second number is the number of words.

#Finds the high frequency words. Can be used to further filter the dtm to remove
##high frequency bias
high_freq_words <- findFreqTerms(tdm, 5000)
high_freq_words
write.table(high_freq_words, file = "highfreq_in_corpus.csv", sep = ";",qmethod = "double")

#Takes the top 40 words and gives counts of each word
tdm_topN <- head(sort(col_sums(tdm), decreasing = TRUE), n=40L)
tdm_topN
write.table(tdm_topN, file = "topN_words_in_corpus.csv", sep = ";",qmethod = "double")

##This takes the term document matrix from the frequncy of word occurance
##and preforms a hierachal cluster analysis on the words

topN.df <- as.data.frame(tdm_topN) #Convert to a data frame 
nrow(topN.df) #number of row
ncol(topN.df) #number of columns

topN.df.scale <- scale(topN.df) #Scales to normalize the data
topN_d <- dist(topN.df.scale, method = "euclidean")
topN_fit <- hclust(topN_d, method="ward")

#Fancy phylogenetic tree version
# pal = needs to have the number of colors equal to the number of your clusters

#Pseduo-ggplot2 version of the denogram
ggdendrogram(topN_fit, theme_dendro=FALSE)
ggdendrogram(topN_fit, rotate=TRUE, size=4, theme_dendro=FALSE, color="tomato") + 
  aes(x=Words, y=Distance) + 
  labs(title="Clustering of the 40 Most Frequent Words") +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,panel.background = element_blank()
  )

pal = c("#556270", "#4ECDC4", "#1B676B", "#FF6B6B", "#C44D58","999999","#C45D58")
#clus takes the number of clusters as the argument fit, #of_clusters
clus = cutree(topN_fit, 7)
op = par(bg="white") #Sets the background color
plot(as.phylo(topN_fit), type="fan",tip.color=pal[clus] ,label.offset=0.3, cex=0.3, col="red")

##Should show something that approximates number of topics present in a corpus
##NOTE: This first part is NOT latent semantic analysis or LDA
td.mat <- as.matrix(tdm)
dist.mat <- dist(t(as.matrix(td.mat)))
dist.mat
fit <- cmdscale(dist.mat, eig = TRUE, k = 2)
points <- data.frame(x = fit$points[, 1], y = fit$points[, 2])
points
ggplot(points, aes(x = x, y = y)) + geom_point(data = points, aes(x = x, y = y)) + 
  labs(title="Multidimensional Scaling of DBR Corpus") +
  aes(x=Coordinate_1, y=Coordinate_2)

k <- 7
SEED <- 2010
jss_TM <-
  list(VEM = LDA(tdm, k = k, control = list(seed = SEED)),
       VEM_fixed = LDA(tdm, k = k,
                       control = list(estimate.alpha = FALSE, seed = SEED)),
       Gibbs = LDA(tdm, k = k, method = "Gibbs",
                   control = list(seed = SEED, burnin = 1000,
                                  thin = 100, iter = 1000)),
       CTM = CTM(tdm, k = k,
                 control = list(seed = SEED,
                                var = list(tol = 10^-4), em = list(tol = 10^-3))))

sapply(jss_TM[1:2], slot, "alpha")

Topic <- topics(jss_TM[["VEM"]], 1)
Terms <- terms(jss_TM[["VEM"]], 10)

topic_list <- Terms[,1:7]
topic_list
write.table(topic_list, file = "topics_in_corpus.csv", sep = ";",qmethod = "double")
topic_list[,1]

(my_topics = topics(jss_TM[["VEM"]])); #Assigns topics to individual documents
topics.df <- as.data.frame(my_topics)
topics.df

ggplot(topics.df, aes(x=my_topics)) + geom_histogram(binwidth=0.5) +
  ggtitle("Distrubtion of Topics in the Corpus") +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,panel.background = element_blank()
  )

#Plot distubtion of topics per document
(d <- Dictionary(c(topic_list[,1:7])))
dict_corp <- inspect(DocumentTermMatrix(rpg_corpus, list(dictionary = d)))
dict.df <- as.data.frame(dict_corp)
dict.df[1:2]
dict.cor <- cor(dict.df,method=c("pearson"))
dict.cor

k4 <- PCA(dict.cor,scale.unit=TRUE, ncp=5, graph=T)
scores = as.data.frame(k4$ind$coord)
colnames(scores) <- c("PCA1", "PCA2")

ggplot(data=scores, aes(x=PCA1, y=PCA2, label=rownames(scores))) + geom_jitter()+
  geom_hline(yintercept=0, colour="gray65") +
  geom_vline(xintercept=0, colour="gray65") +
  geom_text(colour="tomato", alpha=0.8, size=4)+
  ggtitle("PCA of Topic Models") +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,panel.background = element_blank()
  )

# Visualize topic correlations using a force-directed network graph
# network diagram using Fruchterman & Reingold algorithm
library(igraph)
g <- as.undirected(graph.adjacency(dict.cor, weighted=TRUE, mode = "upper"))
# define layout
layout1 <- layout.fruchterman.reingold(g, niter=500)
# use degree to control vertex and font size
b1 <- degree(g) 
V(g)$label.cex <-  b1 * 0.5 / max(b1) # label text size
V(g)$size <-  b1 * 15 / max(b1)        # node size
plot(g, layout = layout1, edge.curved = TRUE, 
     vertex.color = "white", 
     edge.arrow.size = 0, 
     vertex.label.family = 'sans'
)
