#
# DSEA: the Second Step - ENRICHMENT
#
# Author: Dmitrii Bychkov, FIMM 2013
# (dmitrii.bychkov@helsinki.fi)
#######################################################

setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('pipeline_sup.R')

library(grid)
library(gplots)
library(ggplot2)
library(reshape2)
require(Nozzle.R1)

# 1. Upload a New Screen
read.csv(file="../datasets/all_leukemia_cl_june_13_disha_astrid.csv", head=TRUE, sep=",") -> leukemia.DATA

# 2. Identify (a) top Sensitive and (b) Resistant Drugs
cell.line <- "SR"

plotDrugSensitivity(leukemia.matrix, cell.line)
drugs.sensitive <- topSensitive(leukemia.matrix, cell.line)

# 3. Upload corresponding data set wit clusters
load('leukemiaClust.RData')

# 4. Push Both Sets for Enrichment
# That is to verify that most of sensitive drugs 
# from a set tend to appear in the same cluster

enrichment.table <- data.frame(Cluster=character(), Cluster.Size=numeric(), DB.Size=numeric(), Sensitive.Count=numeric(), Sensitive.Total=numeric(), p.Value=numeric(), stringsAsFactors=FALSE)
for(cluster in unique(leukemia.ClUST$Cluster)){
  cluster.set <- leukemia.ClUST[ leukemia.ClUST[,"Cluster"] == cluster, "DrugName" ]
  
  enrichment.table[cluster,"Cluster"] <- paste("Cluster ", cluster, sep="")
  enrichment.table[cluster,"Cluster.Size"] <- length(cluster.set)
  enrichment.table[cluster,"DB.Size"] <- length( unique(leukemia.ClUST$DrugName) ) 
  enrichment.table[cluster,"Sensitive.Count"] <- sum(drugs.sensitive$DrugName %in% cluster.set)
}

