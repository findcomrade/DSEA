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
read.csv(file="../datasets/all_leukemia_cl_june_13_disha_astrid.csv", head=TRUE, sep=",") -> data.NEW

# 2. Identify (a) top Sensitive and (b) Resistant Drugs
target.cell.line <- "THP.1"
threshold <- 0.05            # higher p-value

matrix.NEW <- data.NEW
matrix.NEW <- data.matrix(matrix.NEW[,-c(1,2)])  # del 1st & 2nd rows 
rownames(matrix.NEW) <- data.NEW[,2]  # assign colnames with drug names
drop <- which(apply(matrix.NEW,1,sum) == 0)
matrix.NEW <- matrix.NEW[-drop,]
nas <- is.na(matrix.NEW); matrix.NEW[nas] <- 0
remove(nas, drop)

drugSensitivity(matrix.NEW, target.cell.line)
plot( density( matrix.NEW[,target.cell.line], na.rm=TRUE), main = "Full Set", xlab = "DSS" )
hist(matrix.NEW[,target.cell.line])

drugs.sensitive <- topSensitive(matrix.NEW, target.cell.line, 19)
drugs.resistant <- topResistant(matrix.NEW, target.cell.line)

target.top.count <- length(drugs.sensitive)

# 3. Upload corresponding data set wit clusters
load('RData/leukemiaClust.RData')


# 4. Claculate enrichment
# target.CL - the one we paste to get enrichment
# reference.CL - all those that are already available in a data set 
enrichment.table <- buildEnrichmentCL(matrix.AML, matrix.NEW, target.cell.line)  # (leukemia.ClUST, drugs.sensitive, drugs.resistant) 

dss.sample <- matrix.NEW[,target.cell.line]
correlation.table <- correlationTable(matrix.AML, dss.sample, target.cell.line)

n <- 6
drop <- which( correlation.table[1:n,"Cell.Line"] %in% enrichment.table[1:n,"Cell.Line"] )
tmp <- as.matrix(correlation.table[drop,])


