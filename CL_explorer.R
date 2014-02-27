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

ggSet <- matrix.MERGED #[210:240,]
mggSet <- melt(ggSet)
colnames(mggSet) <- c("Drug", "Sample", "DSS")
ggplot(mggSet, aes(reorder(Drug, DSS, FUN=median, na.rm = TRUE), DSS, fill=Drug)) +  xlab("Drug") +
                      geom_boxplot() + opts(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggSet <- matrix.MERGED #[,1:35]
mggSet <- melt(ggSet)
colnames(mggSet) <- c("Drug", "Sample", "DSS")
ggplot(mggSet, aes(reorder(Sample, DSS, FUN=median, na.rm = TRUE), DSS, fill=Sample)) + xlab("Sample") +
                      geom_boxplot() + opts(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

## CELL LINES
mean.vect <- c()
for(item in colnames(matrix.MERGED)){
  mean.vect <- c( mean.vect, mean(matrix.MERGED[,item], na.rm=TRUE) ) 
}
meanframe <- data.frame(cell.line=colnames(matrix.MERGED), mean.val=mean.vect)
meanframe <- meanframe[ order( meanframe[,"mean.val"], meanframe[,"cell.line"] ), ]

pdf(file='/home/comrade/Desktop/plot.pdf')
boxplot(matrix.MERGED, col=1:length(matrix.MERGED), xlab = "", ylab = "DSS Value", names=NA, at=rank(mean.vect) )
text(seq(1, length(colnames(matrix.MERGED)), by=1), -6, srt = 90, pos = 1, xpd = TRUE, labels=meanframe$cell.line, cex=0.59)
dev.off()

## DRUGS
transponsed.MERGED <- transponsed.MERGED[,210:249]
mean.vect <- c()
for(item in colnames(transponsed.MERGED)){
  mean.vect <- c( mean.vect, mean(transponsed.MERGED[,item], na.rm=TRUE) ) 
}
meanframe <- data.frame(cell.line=colnames(transponsed.MERGED), mean.val=mean.vect)
meanframe <- meanframe[ order( meanframe[,"mean.val"], meanframe[,"cell.line"] ), ]

pdf(file='/home/comrade/Desktop/plot.pdf')
boxplot(transponsed.MERGED, col=1:length(transponsed.MERGED), xlab = "", ylab = "DSS Value", names=NA, at=rank(mean.vect) )
text(seq(1, length(colnames(transponsed.MERGED)), by=1), -5, srt = 90, pos = 1, xpd = TRUE, labels=meanframe$cell.line, cex=0.59)
dev.off()

## By Cell Line
drugSensitivity(matrix.MERGED, cell.line, 25)
plot( density( matrix.MERGED[,cell.line], na.rm=TRUE), main = paste(cell.line, " DSS Distribution", sep=""), xlab = "DSS" )

# Full Set Density
plot( density( transponsed.MERGED, na.rm=TRUE), main = "Full Set", xlab = "DSS" )

drugs.sensitive <- topSensitive(matrix.MERGED, cell.line, 19)
