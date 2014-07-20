setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('pipeline_sup.R')
source('source/dsea_aux.R')

library(vcd)
library(grid)
library(Hmisc)
library(gplots)
library(stringr)
library(graphics)
library(lattice)
library(reshape2)
library(RJSONIO)

### Input Parameters                 ###
# ==================================== #
dsrt.dataset.file <- "../datasets/merged_dss_new.csv"
sample.to.exclude <- "MOLM.13" 

dss.min  <- 3             # minimal dss value which is appropriate for further analysis
dss.min.fract  <- 0.9     #

draw.heat  <- 0

### Import Data Sets                 ###
# ==================================== #
read.csv(dsrt.dataset.file, head=TRUE, sep="\t") -> dsrt.DATA

#Filter out drugs screened over a few cell lines
filt <- apply(as.matrix(dsrt.DATA[,-c(1,2)]), 1, function(x) sum(is.na(x)) < 0.70 * length(x))  #sum(x < dss.min) < dss.min.fract * length(x))
dsrt.DATA  <- dsrt.DATA[filt,]

# Exclude a Sample
tmp <- which( as.character(colnames(dsrt.DATA)) == sample.to.exclude ) 
dsrt.DATA <- dsrt.DATA[,-tmp]

#remove(dsrt.dataset.file, tmp)
### Data Set Preprocessing           ###
# ==================================== #

# Convert the data to 'double matrix' and: 
matrix.CSamples <- dsrt.DATA  # copy
matrix.CSamples <- data.matrix(matrix.CSamples[,-c(1,2)])      # del 1st & 2nd rows 
matrix.CDrugs <- t(matrix.CSamples)                            # and transform
colnames(matrix.CDrugs) <- dsrt.DATA[,2]                       # assign colnames with drug names
rownames(matrix.CSamples) <- dsrt.DATA[,2]                     # assign colnames with drug names

# Clean (trim) batch.ids 
dsrt.DATA$FIMM.Batch.ID <- dsrt.DATA[,"ID.Drug"]
dsrt.DATA[,"ID.Drug"] <- strtrim(dsrt.DATA[,"ID.Drug"],10)


### Compute Correlation              ###
# ==================================== #

# < Samples >
corr.Sampl <- rcorr(matrix.CSamples)
corr.Sampl <- corr.Sampl$r

# < Drugs >
corr.Drugs <- rcorr(matrix.CDrugs)
corr.Drugs <- corr.Drugs$r
nas <- is.na(corr.Drugs); corr.Drugs[nas] <- 0


### Hierarchial Clustering           ###
# ==================================== #
# Alailable methods for HCLUST: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
drugs.Dist  <- as.dist(1-corr.Drugs)
dendro.cut  <- 2
hr          <- hclust(drugs.Dist, method="ward")
mycl        <- cutree(hr, h=dendro.cut)  #h=max(hr$height)/1.72) 
set.count   <- length(unique(mycl))

auxDendroPlot(hr, mycl, set.count, dendro.cut)
# Compare two clustering results
table(mycl, cutree(hclust(as.dist(1-corr.Drugs), method="ward"), k=15, h=3) )

sampl.Dist  <- as.dist(1-abs(corr.Sampl))
hc <- hclust(sampl.Dist, method="complete")
mysl <- cutree(hc, h=0.5)  #h=max(hc$height)/1.05) 



tree.Drugs <- as.data.frame(mycl)         # to use 'tree.Drugs' in further analysis
tree.Sampl <- as.data.frame(mysl)         # to use 'tree.Sampl' in further analysis

mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
mycolhs <- rainbow(length(unique(mysl)), start=0.1, end=0.9); mycolhs <- mycolhs[as.vector(mysl)] 

myheatcol <- topo.colors(75)


if (draw.heat){
  heatmap.2(matrix.CDrugs, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hr), col=myheatcol, 
            scale="row", trace="none", RowSideColors=mycolhs, ColSideColors=mycolhc, cexRow=0.3, cexCol=0.5) 
  dev.off()  
}

remove(myheatcol, mycolhc, mycolhs)

tree.Drugs  <- preprocessClusters(tree.Drugs, space="drugs")
tree.Sampl  <- preprocessClusters(tree.Sampl, space="samples")

drug.partition.list  <- getClusterContent(tree.Drugs, space="drugs")
smpl.partition.list  <- getClusterContent(tree.Sampl, space="samples")

# Save R Objects to a file
save(tree.Drugs, tree.Sampl, dsrt.DATA, matrix.CSamples, matrix.CDrugs, drug.partition.list, smpl.partition.list, file = "RData/Clusters.RData")

### Save JSON for D3 Clusters       ###
# ==================================== #
tree.Drugs[,"isTop"] <- 0 
dropJSON(tree.Drugs, path='Results/json/drug_clust.json')
tree.Sampl[,"isTop"] <- 0 
colnames(tree.Sampl) <- c("Cluster", "DrugName","isTop")
dropJSON(tree.Sampl, path='Results/json/sample_clust.json')


### Save JSON for D3 Dendogram       ###
# ==================================== #
dendo.list <- hclustToTree(hr, mycl)
jsonTree <- toJSON(dendo.list)
sink('Results/json/drug_dendo.json')
cat( substr(jsonTree, 4, nchar(jsonTree)-3) )
sink()

dendo.list <- hclustToTree(hc, mysl)
jsonTree <- toJSON(dendo.list)
sink('Results/json/sample_dendo.json')
cat( substr(jsonTree, 4, nchar(jsonTree)-3) )
sink()

remove(dendo.list, jsonTree)
# ==================================== #

