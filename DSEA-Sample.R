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
library(Hmisc)

source('source/dsea_aux.R')
### Input Parameters                 ###
# ==================================== #
dsrt.dataset.file     <- "../datasets/merged_dss_new.csv"
target.cell.line      <- "SR"
p.cutoff              <- 0.05
cor.cutoff            <- 0.75
dss.cutoff            <- 21
# ==================================== #


# 1. Upload a New Screen
read.csv(dsrt.dataset.file, head=TRUE, sep="\t") -> dsrt.DATA

matrix.New <- dsrt.DATA
matrix.New <- data.matrix(matrix.New[,-c(1,2)])     # del 1st & 2nd rows 
rownames(matrix.New) <- dsrt.DATA[,2]               # assign colnames with drug names

# ==================================== #

drugSensitivity(matrix.New, target.cell.line, dss.cutoff)

sample.profile <- matrix.New['Imatinib',]
prob.df <- auxDSSDensityEstimate(vector=sample.profile, hh.cells=55, graph=TRUE,
                          title=paste("Sample", target.cell.line, "PDF Estimate", sep=" :: "), xax.text="DSS")

auxDSSSpecificityScore(prob.df, 5, graph=TRUE)

drugs.sensitive <- topSensitive(matrix.New, target.cell.line, dss.cutoff)
drugs.resistant <- topResistant(matrix.New, target.cell.line)

target.top.count <- length(drugs.sensitive$DrugName)


# 3. Upload corresponding data set wit clusters
load('RData/Clusters.RData')

target.list    <- drugRankedList(matrix.New, target.cell.line)
reference.list <- drugRankedList(matrix.CSamples, "SIG.M5")


# 4. Claculate enrichment
# target.CL - the one we paste to get enrichment
# reference.CL - all those that are already available in a data set 
enrichment.table <- buildEnrichmentCL(matrix.CSamples, matrix.New, target.cell.line, dss.cutoff)

dss.sample <- matrix.New[,target.cell.line]
correlation.table <- correlationTable(matrix.CSamples, dss.sample, target.cell.line)

enriched.cell.lines <- getEnrichedCellLines(correlation.table, enrichment.table, cor.cutoff, p.cutoff)

# Add information (new col) to 'tree.Samples' about
# which sample to highlight: as enriched
is.top <- tree.Sampl[,"SampleName"] %in% enriched.cell.lines[,"Cell.Line"]

tree.Sampl[,"isTop"] <- 0 
tree.Sampl[is.top,"isTop"] <- 1
colnames(tree.Sampl) <- c("Cluster", "DrugName","isTop")

dropJSON(tree.Sampl, path='Results/json/sample_clust.json')
