setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('pipeline_sup.R')

library(grid)
library(gplots)
library(ggplot2)
library(reshape2)

load('leukemiaClust.RData')
# list cell.lines from dataset
colnames(leukemia.matrix)
cell.line <- "SIG.M5"

## CELL LINES
mean.vect <- c()
for(item in colnames(leukemia.matrix)){
  mean.vect <- c( mean.vect, mean(leukemia.matrix[,item]) ) 
}
meanframe <- data.frame(cell.line=colnames(leukemia.matrix), mean.val=mean.vect)
meanframe <- meanframe[ order( meanframe[,"mean.val"], meanframe[,"cell.line"] ), ]

boxplot(leukemia.matrix, col=1:length(leukemia.matrix), xlab = "Cell Lines", ylab = "DSS Value Distribution", names=NA, at=rank(mean.vect) )
text(seq(1, length(colnames(leukemia.matrix)), by=1), -4, srt = 90, pos = 1, xpd = TRUE, labels=meanframe$cell.line, cex=0.59)

## DRUGS
mean.vect <- c()
for(item in colnames(leukemia.transponsed)){
  mean.vect <- c( mean.vect, mean(leukemia.transponsed[,item]) ) 
}
meanframe <- data.frame(cell.line=colnames(leukemia.transponsed), mean.val=mean.vect)
meanframe <- meanframe[ order( meanframe[,"mean.val"], meanframe[,"cell.line"] ), ]

boxplot(leukemia.transponsed, col=1:length(leukemia.transponsed), xlab = "Cell Lines", ylab = "DSS Value Distribution", names=NA, at=rank(mean.vect) )
text(seq(1, length(colnames(leukemia.transponsed)), by=1), -4, srt = 90, pos = 1, xpd = TRUE, labels=meanframe$cell.line, cex=0.59)

## By Cell Line
drugSensitivity(leukemia.matrix, cell.line, 25)
plot( density( leukemia.matrix[,cell.line], na.rm=TRUE), main = paste(cell.line, " DSS Distribution", sep=""), xlab = "DSS" )


drugs.sensitive <- topSensitive(leukemia.matrix, cell.line, 19)
