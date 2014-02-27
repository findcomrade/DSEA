setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('pipeline_sup.R')

library(vcd)
library(grid)
library(Hmisc)
library(gplots)
library(stringr)
library(lattice)
library(reshape2)
library(RJSONIO)

### Input Parameters                 ###
# ==================================== #
dsrt.dataset.file <- "../datasets/merged_dss_new.csv"
sample.to.exclude <- "SR"

### Import Data Sets                 ###
# ==================================== #
read.csv(dsrt.dataset.file, head=TRUE, sep="\t") -> dsrt.DATA

# Exclude a Sample
tmp <- which( colnames(dsrt.DATA) == sample.to.exclude ) 
dsrt.DATA <- dsrt.DATA[,-tmp]

remove(dsrt.dataset.file, tmp)
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
hr <- hclust(as.dist(1-corr.Drugs), method="complete")
hc <- hclust(as.dist(1-corr.Sampl), method="complete")  

# Cuts the tree and creates color vector for clusters
mycl <- cutree(hr, h=max(hr$height)/1.5) 
mysl <- cutree(hc, h=max(hr$height)/4) 

tree.Drugs <- as.data.frame(mycl)         # to use 'tree.Drugs' in further analysis
tree.Sampl <- as.data.frame(mysl)         # to use 'tree.Sampl' in further analysis

mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
mycolhs <- rainbow(length(unique(mysl)), start=0.1, end=0.9); mycolhs <- mycolhs[as.vector(mysl)] 

myheatcol <- topo.colors(75)

# Creates heatmap for entire data set where the obtained clusters are indicated in the color bar
heatmap.2(matrix.CSamples, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, 
          scale="row", trace="none", RowSideColors=mycolhc, ColSideColors=mycolhs, cexRow=0.3, cexCol=0.5) 
dev.off()

remove(myheatcol, mycolhc, mycolhs)

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


# Prepare Obtained Clusters for Drugs
tree.Drugs$DrugName <- rownames(tree.Drugs)
colnames(tree.Drugs) <- c("Cluster", "DrugName")
tree.Drugs$DrugName <- lapply(tree.Drugs[,"DrugName"], str_trim)

# Prepare Obtained Clusters for Samples
tree.Sampl$SampleName <- rownames(tree.Sampl)
colnames(tree.Sampl) <- c("Cluster", "SampleName")
tree.Sampl$SampleName <- lapply(tree.Sampl[,"SampleName"], str_trim)

# Get Drug Lists by Clusters
drug.partition.list <- data.frame()
for( clst in unique(tree.Drugs$Cluster) ){
  drug.partition.list[clst,1] <- paste("Cluster ", clst, sep="")
  drug.partition.list[clst,2] <- length(tree.Drugs[tree.Drugs[,"Cluster"] == clst,"DrugName"])
  drug.partition.list[clst,3] <- gsub(",", "", toString( tree.Drugs[tree.Drugs[,"Cluster"] == clst,"DrugName"] )) 
}
colnames(drug.partition.list) <- c("ClusterSym", "Cluster.Size", "Drug.List")

# Get Sample Lists by Clusters
smpl.partition.list <- data.frame()
for( clst in unique(tree.Sampl$Cluster) ){
  smpl.partition.list[clst,1] <- paste("Cluster ", clst, sep="")
  smpl.partition.list[clst,2] <- length(tree.Sampl[tree.Sampl[,"Cluster"] == clst,"SampleName"])
  smpl.partition.list[clst,3] <- gsub(",", "", toString( tree.Sampl[tree.Sampl[,"Cluster"] == clst,"SampleName"] )) 
}
colnames(smpl.partition.list) <- c("ClusterSym", "Cluster.Size", "Sample.List")

# Save R Objects to a file
save(tree.Drugs, tree.Sampl, dsrt.DATA, matrix.CSamples, matrix.CDrugs, drug.partition.list, smpl.partition.list, file = "RData/Clusters.RData")


