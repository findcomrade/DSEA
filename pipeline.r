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
library(RJSONIO)
library(reshape2)
require(Nozzle.R1)

# 1. Upload a New Screen
read.csv(file="../datasets/all_leukemia_cl_june_13_disha_astrid.csv", head=TRUE, sep=",") -> leukemia.DATA

# 2. Identify (a) top Sensitive and (b) Resistant Drugs
cell.line <- "SR"

leukemia.matrix <- leukemia.DATA
leukemia.matrix <- data.matrix(leukemia.matrix[,-c(1,2)])  # del 1st & 2nd rows 
rownames(leukemia.matrix) <- leukemia.DATA[,2]  # assign colnames with drug names
drop <- which(apply(leukemia.matrix,1,sum) == 0)
leukemia.matrix <- leukemia.matrix[-drop,]
nas <- is.na(leukemia.matrix); leukemia.matrix[nas] <- 0

drugSensitivity(leukemia.matrix, cell.line)
plot( density( leukemia.matrix[,cell.line], na.rm=TRUE), main = "Full Set", xlab = "DSS" )
hist(leukemia.matrix[,cell.line])

drugs.sensitive <- topSensitive(leukemia.matrix, cell.line, 19)
drugs.resistant <- topResistant(leukemia.matrix, cell.line)

# 3. Upload corresponding data set wit clusters
load('leukemiaClust.RData')

# 4. Push Both Sets for Enrichment
# That is to verify that most of sensitive drugs 
# from a set tend to appear in the same cluster

enrichment.table <- buildEnrichment(leukemia.ClUST, drugs.sensitive, drugs.resistant)

# Add information (new col) to 'leukemia.ClUST' about
# which drugs to highlight: sensitive or resistant
is.top <- leukemia.ClUST[,"DrugName"] %in% drugs.sensitive$DrugName # sensit
is.bot <- leukemia.ClUST[,"DrugName"] %in% drugs.resistant$DrugName # resist
leukemia.ClUST[,"isTop"] <- 0 
leukemia.ClUST[is.top,"isTop"] <- 1
leukemia.ClUST[is.bot,"isTop"] <- -1

dropJSON(leukemia.ClUST)

