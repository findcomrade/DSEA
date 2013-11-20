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

# 1. Upload a New Screen
read.csv(file="../datasets/merged_dss_new.csv", head=TRUE, sep="\t") -> data.NEW

# 2. Identify (a) top Sensitive and (b) Resistant Drugs
target.cell.line <- "primary_culture_270513_PCA2_ROCK"
threshold <- 0.05            # higher p-value

matrix.NEW <- data.NEW
matrix.NEW <- data.matrix(matrix.NEW[,-c(1,2)])  # del 1st & 2nd rows 
rownames(matrix.NEW) <- data.NEW[,2]  # assign colnames with drug names
#drop <- which(apply(matrix.NEW,1,sum) == 0)
#matrix.NEW <- matrix.NEW[-drop,]
#nas <- is.na(matrix.NEW); matrix.NEW[nas] <- 0
remove(nas, drop)

drugSensitivity(matrix.NEW, target.cell.line, 15)
# plot( density( matrix.NEW[,target.cell.line], na.rm=TRUE), main = "Full Set", xlab = "DSS" )
# hist(matrix.NEW[,target.cell.line])

dss.cutoff <- 9
drugs.sensitive <- topSensitive(matrix.NEW, target.cell.line, dss.cutoff)
drugs.resistant <- topResistant(matrix.NEW, target.cell.line)

target.top.count <- length(drugs.sensitive)

# 3. Upload corresponding data set wit clusters
load('RData/leukemiaClust.RData')


# 4. Claculate enrichment
# target.CL - the one we paste to get enrichment
# reference.CL - all those that are already available in a data set 
enrichment.table <- buildEnrichmentCL(matrix.MERGED, matrix.NEW, target.cell.line, dss.cutoff)  # (leukemia.ClUST, drugs.sensitive, drugs.resistant) 

dss.sample <- matrix.NEW[,target.cell.line]
correlation.table <- correlationTable(matrix.MERGED, dss.sample, target.cell.line)

enriched.cell.lines <- getEnrichedCellLines(correlation.table, enrichment.table, 0.75, 0.05)

# Add information (new col) to 'tree.DRUGS' about
# which drugs to highlight: sensitive or resistant
is.top <- tree.SAMPLES[,"SampleName"] %in% enriched.cell.lines[,"Cell.Line"]
# enrichment.table[1:6, "Cell.Line"] # correlation.table[1:21,"Cell.Line"] # enriched.cell.lines$Cell.Line  

tree.SAMPLES[,"isTop"] <- 0 
tree.SAMPLES[is.top,"isTop"] <- 1
colnames(tree.SAMPLES) <- c("Cluster", "DrugName","isTop")

dropJSON(tree.SAMPLES, path='/home/comrade/Projects/d3.v3/CL-tree.json')
