require(ade4)
require(stats)
library(Hmisc)
require(ggplot2)
require(reshape2)
require(FactoMineR)

getDrugsSensitivityRank <- function(cDrug.matrix){
  #  Produce Sensitivity Profile for each drug in the matrix 
  #  and order drugs according to perc of high sensitivity.
  #
  #  cDrug.matrix:  double matrix with Drugs as rownames;
  
  xpr.melt <- melt(cDrug.matrix)
  fct <- cut(xpr.melt$value, breaks = c(min(xpr.melt$value,na.rm=TRUE) * 2,5,15, max(xpr.melt$value,na.rm=TRUE) * 2))
  levels(fct) <- c("Resistant (< 5)","Intermediate (5 - 15)","Sensitive (>= 15)")
  mat <-  table(xpr.melt$value,fct)
  partition <-  prop.table(mat,1)
  partition <- partition[,3:1]
  partition <- partition[order(partition[,1],decreasing=TRUE),]
  #partition$Rank <- c(1:nrow(partition))
    
  return(partition)
}

plotPCAScores <- function(scores, header){
  # plotPACScores(pca.prcomp$x, "PCA")
  scores = as.data.frame(scores)
  
  ggplot(data=scores, aes(x=PC1, y=PC2, label=rownames(scores))) +
                        geom_hline(yintercept=0, colour="gray65") +
                        geom_vline(xintercept=0, colour="gray65") +
                        geom_text(colour="tomato", alpha=0.8, size=4) +
                        opts(title=header)
}

getRank <- function(drug.name){
  return( which(rownames(partition) == drug.name) )
}