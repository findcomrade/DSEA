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

### Input Parameters                 ###
# ==================================== #
dsrt.dataset.file     <- "../datasets/merged_dss_new.csv"
target.sample         <- "SR"
dss.cutoff            <- 21
# ==================================== #


# 1. Upload a New Screen
read.csv(dsrt.dataset.file, head=TRUE, sep="\t") -> dsrt.DATA

matrix.CSamples <- dsrt.DATA
matrix.CSamples <- data.matrix(matrix.CSamples[,-c(1,2)])    # del 1st & 2nd rows 
rownames(matrix.CSamples) <- dsrt.DATA[,2]                   # assign colnames with drug names

# ==================================== #

drugSensitivity(matrix.CSamples, target.sample, dss.cutoff)
sample.profile <- matrix.CSamples[,target.sample]

drugs.sensitive <- topSensitive(matrix.CSamples, target.sample, dss.cutoff)
drugs.resistant <- topResistant(matrix.CSamples, target.sample)

drug.list <- drugs.sensitive$DrugName
x <- matrix.CSamples[drug.list,]

# 3. Upload corresponding data set wit clusters
load('RData/Clusters.RData')


# 4. Push Both Sets for Enrichment
# That is to verify that most of sensitive drugs 
# from a set tend to appear in the same cluster
enrichment.table <- buildEnrichmentD(tree.Drugs, drugs.sensitive, drugs.resistant)

# Add information (new col) to 'tree.DRUGS' about
# which drugs to highlight: sensitive or resistant
is.top <- tree.Drugs[,"DrugName"] %in% drugs.sensitive$DrugName      # sensitive
is.bot <- tree.Drugs[,"DrugName"] %in% drugs.resistant$DrugName      # resistant
tree.Drugs[,"isTop"] <- 0 
tree.Drugs[is.top,"isTop"] <- 1
#tree.DRUGS[is.bot,"isTop"] <- -1

dropJSON(tree.Drugs, path='Results/json/drug_clust.json')
