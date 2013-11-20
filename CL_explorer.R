#
# This is cell line explorer.
# You can plot DSS values distrubution
# across: a) Cell Lines b) Drugs
# as well as estimate sensitivity behaviour 
# of a cell line by a bar plot for druds DSS values
#
#######################################################

setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('pipeline_sup.R')

library(grid)
library(gplots)
library(ggplot2)
library(reshape2)

load('RData/leukemiaClust.RData')
# list cell.lines from dataset
colnames(matrix.MERGED)
cell.line <- "SIG.M5"

## CELL LINES
mean.vect <- c()
for(item in colnames(matrix.MERGED)){
  mean.vect <- c( mean.vect, mean(matrix.MERGED[,item]) ) 
}
meanframe <- data.frame(cell.line=colnames(matrix.MERGED), mean.val=mean.vect)
meanframe <- meanframe[ order( meanframe[,"mean.val"], meanframe[,"cell.line"] ), ]

boxplot(matrix.MERGED, col=1:length(matrix.MERGED), xlab = "Cell Lines", ylab = "DSS Value Distribution", names=NA, at=rank(mean.vect) )
text(seq(1, length(colnames(matrix.MERGED)), by=1), -4, srt = 90, pos = 1, xpd = TRUE, labels=meanframe$cell.line, cex=0.59)

## DRUGS
mean.vect <- c()
for(item in colnames(transponsed.MERGED)){
  mean.vect <- c( mean.vect, mean(transponsed.MERGED[,item]) ) 
}
meanframe <- data.frame(cell.line=colnames(transponsed.MERGED), mean.val=mean.vect)
meanframe <- meanframe[ order( meanframe[,"mean.val"], meanframe[,"cell.line"] ), ]

boxplot(transponsed.MERGED, col=1:length(transponsed.MERGED), xlab = "Cell Lines", ylab = "DSS Value Distribution", names=NA, at=rank(mean.vect) )
text(seq(1, length(colnames(transponsed.MERGED)), by=1), -4, srt = 90, pos = 1, xpd = TRUE, labels=meanframe$cell.line, cex=0.59)

## By Cell Line
drugSensitivity(matrix.MERGED, cell.line, 25)
plot( density( matrix.MERGED[,cell.line], na.rm=TRUE), main = paste(cell.line, " DSS Distribution", sep=""), xlab = "DSS" )

# Full Set Density
plot( density( transponsed.MERGED, na.rm=TRUE), main = "Full Set", xlab = "DSS" )

drugs.sensitive <- topSensitive(matrix.MERGED, cell.line, 19)
