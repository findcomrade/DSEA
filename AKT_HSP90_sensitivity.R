# Aquired sensitivity for PI3K/AKT and HSP90 inhibitors
# in resistant variants in comparison to parental
#
# Author: Dmitrii Bychkov, FIMM 2013
#######################################################

setwd("/home/comrade/Ubuntu One/DSEA/r-code")


library(grid)
library(gplots)
library(ggplot2)
library(RJSONIO)
library(reshape2)

### Input Parameters                 ###
# ==================================== #
dsrt.dataset.file     <- "../datasets/DSS2_merged_34-samples_FO3D-16-onwards_2014-04-15.csv"
target.sample         <- "LNCaP.AF8.210113"
dss.cutoff            <- 11
# ==================================== #


# 1. Upload a New Screen
read.csv(dsrt.dataset.file, head=TRUE, sep="\t") -> dsrt.DATA
dsrt.DATA  <- dsrt.DATA[,c(2,23:32)]

drug.lst  <- c("Tanespimycin", "XL147", "MK+AC0-2206", "AZD8055", "GDC+AC0-0980", "XL765", "Geldanamycin")

# Filter drugs:
ind  <- which( dsrt.DATA[,"Name.Drug"] %in% drug.lst )
dsrt.DATA  <- dsrt.DATA[ind,]

# PC346C Flu1
tbl  <- dsrt.DATA[,c('Name.Drug', 'PC346C.AF8.181113', 'PC346C.AF8.Flu1.AF8.181113')]
# PC346C Flu2
tbl  <- dsrt.DATA[,c('Name.Drug', 'PC346C.AF8.181113', 'PC346.AF8.Flu2.AF8.111113')]
# PC346C DCC
tbl  <- dsrt.DATA[,c('Name.Drug', 'PC346C.AF8.181113', 'PC346.AF8.DCC.AF8.111113')]

# VCaP rep1 BIC
tbl  <- dsrt.DATA[,c('Name.Drug', 'VCaP.AF8.290413.AF8.replicate1', 'VCaP.BIC.AF8.290413.AF8.replicate1')]
# VCaP rep1 FLU
tbl  <- dsrt.DATA[,c('Name.Drug', 'VCaP.AF8.290413.AF8.replicate1', 'VCaP.FLU.AF8.230513.AF8.replicate1')]
# VCaP rep2 BIC
tbl  <- dsrt.DATA[,c('Name.Drug', 'VCaP.AF8.230513.AF8.replicate2', 'VCaP.BIC.AF8.230513.AF8.replicate2')]



################
tbl  <- melt(tbl)

anno.row  <- rbind( data.frame(status=rep("parental", 7)), data.frame(status=rep("resistant", 7)) )
tbl       <- cbind(tbl, anno.row)
colnames(tbl)  <- c("Drug.Name", "Sample", "DSS", "Status")

ggplot(data=tbl, aes(x=factor(Drug.Name), y=DSS, fill=Status)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  ggtitle("VCaP_replicate2\nBIC resistant clone")


