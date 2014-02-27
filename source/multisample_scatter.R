library(Hmisc)
library(ggplot2)
library(reshape2)
library(ggdendro)

# A Scatter Plot with IDs on X and Samples (WSNwt, X202...) on Y
# Below is an example of input file to be loaded in expr.data

# ID       WSNwt	    X202	      X210	        X220	      RK.AA
# CXCL10	7,107022349	7,15347294	5,0120528058	7,125860322	5,23475358
# CCL8	  6,022000793	7,24845859	2,7592660878	6,409497927	4,80852799
#

dataset.file <- "/home/comrade/Projects/Kainov/Illumina_HumanHT_12_221013/Anastasina (diff WSN mutants)/macros_cut3_gene_expr_data_for_plotting.csv"
read.csv(dataset.file, head=TRUE, sep="\t") -> expr.data
remove(dataset.file)

expr.matrix             <- data.matrix(expr.data[,2:6])                             # conver data frame to double matrix
rownames(expr.matrix)   <- expr.data[,1]                                            # give names to matrix rows

corr.Genes <- rcorr( t(expr.matrix) )                                               # Compute correlations
corr.Genes <- corr.Genes$r

hc <- hclust(as.dist(1-corr.Genes), method="complete")                              # clustering

dd.col  <- as.dendrogram(hc)                                                        # get order of genes from hclust
col.ord <- order.dendrogram(dd.col)                                                 # get order of genes from hclust
 
#expr.ordered          <- expr.data[col.ord,]                                       #   Construct new matrix with order 
expr.ordered          <- expr.data[ order(-expr.data[,6]), ]                        #   Construct new matrix with order
expr.ordered$ID       <- with(expr.ordered, factor(ID, levels=ID, ordered=TRUE))
expr.melt             <- melt(expr.ordered, id.vars="ID")
colnames(expr.melt)   <- c("Gene", "Sample", "logFC")

ggplot(expr.melt, aes(x=Gene, y=logFC, color=Sample)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, size=5, hjust = 1))

plot(hc)
