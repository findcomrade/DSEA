setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('source/explaratoryAnalysisSource.R')

read.csv(file="../datasets/merged_dss_new.csv", head=TRUE, sep="\t") -> initial.DSRT
read.csv(file="../datasets/Samples.csv", head=TRUE, sep="\t")        -> samples.Supp

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

drugs.pca  <- PCA(matrix.DSRT[,1:183], ncp=8, scale.unit=TRUE, quanti.sup=c(1), graph=TRUE)

samples.pca  <- PCA(tmatrix.DSRT[,1:230], ncp=5, scale.unit=TRUE, quali.sup=c(1,2), graph=TRUE)
pdf(file = "smaples_pca.pdf", width = 15, height = 10,)
plot.PCA(samples.pca, axes=c(1, 2), choix="ind", habillage=1, cex=0.4)
dev.off()

pdf(file = "smaples_pca_zoom.pdf")
plot.PCA(samples.pca, axes=c(1, 2), choix="ind", habillage=1, cex=0.3, xlim=c(-10,20), ylim=c(-5,7))
dev.off()

samples.clustering <- HCPC(samples.pca)


### Plot Loadings                    ###
# ==================================== #
samples.corr <- rcorr(matrix.DSRT[,2:183])$r
qplot(x=X1, y=X2, data=melt(samples.corr), geom="tile",  fill=value)

pca <- prcomp(na.omit(matrix.DSRT[,2:183]), scale=T)
melted <- cbind( samples.Supp$Sample.Origin, melt(pca$rotation[,1:4]) )
colnames(melted)[1] <- "origin"

ggplot(data=melted) + geom_bar(aes(x=X1, y=value, fill=origin), stat="identity") +  facet_wrap(~X2)

### OLD
#======================
# ver 1
# pca.prcomp    <- prcomp(na.omit(subSet), scale. = TRUE)
# biplot(pca.prcomp)
# ver 2
# pca.dudi      <- dudi.pca(na.omit(tmatrix.DSRT), nf=5, scannf=FALSE)
# biplot(pca.dudi)








