setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('source/explaratoryAnalysisSource.R')
source('source/dsea_aux.R')
library(ggplot2)
library(gplots)
library(reshape2)

read.csv(file="../datasets/merged_dss_new.csv", head=TRUE, sep="\t") -> initial.DSRT
read.csv(file="../datasets/Samples.csv", head=TRUE, sep="\t")        -> samples.Supp
# prostate set
# read.csv(file="../datasets/DSS2_merged_34-samples_FO3D-16-onwards_2014-04-15.csv", head=TRUE, sep="\t") -> initial.DSRT

# ==================================== #
### Data Preprocessing               ###
# ==================================== #
filt <- apply(data.matrix(initial.DSRT[,-c(1,2)]), 1, function(x) sum(is.na(x)) < 0.50 * length(x))
initial.DSRT <- initial.DSRT[filt,]                   # Drop Drugs with NAs more than 50%

matrix.DSRT  <- data.matrix(initial.DSRT[,-c(1,2)])   # del 1st (ID) & 2nd (Alias) cols 
tmatrix.DSRT <- t(matrix.DSRT)                        # and transform
rownames(matrix.DSRT)  <- initial.DSRT[,2]            # assign drug names
colnames(tmatrix.DSRT) <- initial.DSRT[,2]            # assign drug names

remove(filt)

# ==================================== #
###      Correlations                ###
# ==================================== #
samples.dss  <- matrix.DSRT[,1]
samples.rank <- as.matrix( rank(samples.dss) )

plot(samples.dss[samples.rank])

# ==================================== #
### PCA :: Supplementary Data        ###
# ==================================== #
# Add Ranks as the 1st column
partition <- getDrugsSensitivityRank(matrix.DSRT)
matrix.DSRT <- cbind(as.numeric(lapply(rownames(matrix.DSRT), getRank)), matrix.DSRT)
colnames(matrix.DSRT)[1] <- "Rank"

#tmatrix.DSRT <- cbind(samples.Supp[,c(2,3)], tmatrix.DSRT)
       
remove(partition, samples.Supp)
# ==================================== #
### PCA                              ###
# ==================================== #

subSet <- matrix.DSRT #[14:27,1:11]

drugs.pca  <- PCA(matrix.DSRT[,2:182], ncp=8, scale.unit=TRUE, graph=TRUE) # quanti.sup=c(1)
samples.pca  <- PCA(tmatrix.DSRT[,1:228], ncp=5, scale.unit=TRUE, quali.sup=c(1,2), graph=TRUE)
plot.PCA(drugs.pca, axes=c(1, 2), choix="ind", cex=0.4)


pdf(file = "figures/smaples_pca.pdf", width = 15, height = 10,)
plot.PCA(samples.pca, axes=c(1, 2), choix="ind", habillage=1, cex=0.4)
dev.off()

pdf(file = "figures/smaples_pca_zoom.pdf")
plot.PCA(samples.pca, axes=c(1, 2), choix="ind", habillage=1, cex=0.3, xlim=c(-10,20), ylim=c(-5,7))
dev.off()

# samples.clustering <- HCPC(samples.pca)


### Plot Loadings                    ###
# ==================================== #
samples.corr <- rcorr(matrix.DSRT[,2:182])$r
qplot(x=X1, y=X2, data=melt(samples.corr), geom="tile",  fill=value)

pca <- prcomp(na.omit(matrix.DSRT[,1:182]), scale=T)
melted <- cbind( samples.Supp$Sample.Origin, melt(pca$rotation[,1:4]) )
colnames(melted)[1] <- "origin"

ggplot(data=melted) + geom_bar(aes(x=Var1, y=value, fill=origin), stat="identity") +  facet_wrap(~Var2)

# END of Loadings                   
# ==================================== #



#### Plot Drug and Sample Profiles   ###
# ==================================== #
drug.names    <- rownames(matrix.DSRT)
sample.names  <- colnames(matrix.DSRT)

# Drug Profiles
PDFPath = "figures/DrugsDSSProfiles.pdf"
pdf(file=PDFPath)  

par(mfrow = c(3,3))
for(drug in drug.names){  
  profile <- matrix.DSRT[drug,]
  control <- round(mean(profile, na.rm=TRUE), digits=2)
  prob.df <- auxDSSDensityEstimate(vector=profile, graph=FALSE)
  if(class(prob.df) == "data.frame"){
    auxDSSSpecificityScore(probdf=prob.df, thresh=control, graph=TRUE, title=paste(drug,"Profile", sep=" :: "))  
  }  
}
dev.off()

# Sample Profiles
PDFPath = "figures/SampleDSSProfiles.pdf"
pdf(file=PDFPath)  

par(mfrow = c(3,3))
for(sample in sample.names){  
  profile <- matrix.DSRT[,sample]
  control <- round(mean(profile, na.rm=TRUE), digits=2)
  prob.df <- auxDSSDensityEstimate(vector=profile, graph=FALSE)
  if(class(prob.df) == "data.frame"){
    auxDSSSpecificityScore(probdf=prob.df, thresh=control, graph=TRUE, title=paste(sample,"Profile", sep=" :: "))  
  }  
}
dev.off()

# END of Drug and Sample Profiles.
# ==================================== #

#### Tissue Specificity Matrix              ###
# ==================================== #
tissue.means      <- auxTissueDSSMeans(matrix.DSRT, samples.Supp)
write.table(tissue.means, file = "means_data.csv")
# Drug Profiles
drug.p    <- "Plicamycin"
tissue.p  <- "prostate" # "ovarian", "leukemia", "prostate", "breast"

profile <- matrix.DSRT[drug.p,]
control <- round(tissue.means[drug.p,tissue.p], digits=2)
if(!is.na(control)){
  prob.df <- auxDSSDensityEstimate(vector=profile, graph=F)
  print(control)
  if(class(prob.df) == "data.frame"){
    auxDSSSpecificityScore(probdf=prob.df, thresh=control, graph=T, title=paste(drug.p,tissue.p, sep=" in "))  
  }  
}

# END of Tissue Specificity.
# ==================================== #

#### Specificity Matrix              ###
# ==================================== #
individual.specificity  <- auxSpecificityMatrix(matrix.DSRT)
tissue.specificity      <- auxTissueSpecificity(matrix.DSRT, samples.Supp)

tissue.count  <- ncol(tissue.specificity)

heat.data  <- tissue.specificity
filt <- which(  apply(heat.data, 1, function(x) sum(is.na(x)) == tissue.count || sum(x, na.rm=TRUE) == 0 )  )
heat.data  <- heat.data[-c(filt),]

## Prepare for heatmap
heat.data  <- data.frame(cbind(rownames(heat.data), heat.data))
write.table(heat.data, file = "specificity_data.csv")

colnames(heat.data)[1] <- "Drug"
heat.data$leukemia <- as.numeric(as.character(heat.data$leukemia))
heat.data$ovarian <- as.numeric(as.character(heat.data$ovarian))
heat.data$prostate <- as.numeric(as.character(heat.data$prostate))
heat.data$breast <- as.numeric(as.character(heat.data$breast))

heat.ordered          <- heat.data[ order(heat.data[,"prostate"]), ]                        #   Construct new matrix with order
heat.ordered$Drug     <- with(heat.ordered, factor(Drug, levels=Drug, ordered=TRUE))
heat.data.m           <- melt(heat.ordered, id.vars="Drug")

colnames(heat.data.m)  <- c("Drug", "Tissue", "Specificity")

PDFPath = "figures/RowSpecificityHeatMap.pdf"
pdf(file=PDFPath, width=150, height=750)  

# HeatMap with raw values
ggplot(heat.data.m, aes(Tissue, Drug) ) +
  geom_tile(aes(fill = Specificity)) + 
  geom_text(aes(fill = heat.data.m$Specificity,label = round(heat.data.m$Specificity, 0)), size=1.5) +
  scale_fill_gradient(low = "white", high = "red")

dev.off()

### Scaled HeatMap
#heat.data[is.na(heat.data)] <- 0
pdf("figures/ScaledSpecificityHeatMap.pdf", paper="a4", width=8, height=11)
heatmap.2(heat.data,
          main = "Tissue Specificity", # heat map title
          density.info="none",
          trace="none",         # turns off trace lines inside the heat map
          #margins=c(3.3,5),       # widens margins around plot
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",            # turn off column clustering
          scale="column",
          cexRow=0.354,
          cexCol=0.85,
          key=TRUE,
          col=bluered(512),
          na.rm=TRUE ) 

dev.off()

# END of Specificity matrix.
# ==================================== #

