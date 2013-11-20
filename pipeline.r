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

# 1. Upload a New Screen
read.csv(file="../datasets/all_leukemia_cl_june_13_disha_astrid.csv", head=TRUE, sep=",") -> data.AML

# 2. Identify (a) top Sensitive and (b) Resistant Drugs
cell.line <- "SR"

matrix.AML <- data.AML
matrix.AML <- data.matrix(matrix.AML[,-c(1,2)])  # del 1st & 2nd rows 
rownames(matrix.AML) <- data.AML[,2]  # assign colnames with drug names
drop <- which(apply(matrix.AML,1,sum) == 0)
matrix.AML <- matrix.AML[-drop,]
nas <- is.na(matrix.AML); matrix.AML[nas] <- 0
remove(nas, drop)

drugSensitivity(matrix.AML, cell.line)
plot( density( matrix.AML[,cell.line], na.rm=TRUE), main = "Full Set", xlab = "DSS" )
hist(matrix.AML[,cell.line])

drugs.sensitive <- topSensitive(matrix.AML, cell.line, 19)
drugs.resistant <- topResistant(matrix.AML, cell.line)

# 3. Upload corresponding data set wit clusters
load('RData/leukemiaClust.RData')


# 4. Push Both Sets for Enrichment
# That is to verify that most of sensitive drugs 
# from a set tend to appear in the same cluster
enrichment.table <- buildEnrichmentD(leukemia.ClUST, drugs.sensitive, drugs.resistant)

# Add information (new col) to 'leukemia.ClUST' about
# which drugs to highlight: sensitive or resistant
is.top <- leukemia.ClUST[,"DrugName"] %in% drugs.sensitive$DrugName # sensit

is.bot <- leukemia.ClUST[,"DrugName"] %in% drugs.resistant$DrugName # resist
leukemia.ClUST[,"isTop"] <- 0 
leukemia.ClUST[is.top,"isTop"] <- 1
leukemia.ClUST[is.bot,"isTop"] <- -1

dropJSON(leukemia.ClUST, path='/home/comrade/Projects/d3.v3/circular.json')
dropCirclePackingJSON(leukemia.ClUST, path='/home/comrade/Projects/d3.v3/circle_packing.json')
