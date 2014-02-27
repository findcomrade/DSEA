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

# 1. Upload a New Screen
read.csv(file="../datasets/merged_dss_new.csv", head=TRUE, sep="\t") -> data.MERGED  # use check.names
#read.csv(file="../datasets/leukemia_all_celllines_data_DSS.csv", head=TRUE, sep="\t") -> data.MERGED

# 2. Identify (a) top Sensitive and (b) Resistant Drugs
cell.line <- "SR"

matrix.MERGED <- data.MERGED
matrix.MERGED <- data.matrix(matrix.MERGED[,-c(1,2)])  # del 1st & 2nd rows 
rownames(matrix.MERGED) <- data.MERGED[,2]  # assign colnames with drug names
#drop <- which(apply(matrix.MERGED,1,sum) == 0)
#matrix.MERGED <- matrix.MERGED[-drop,]
#nas <- is.na(matrix.MERGED); matrix.MERGED[nas] <- 0
remove(nas, drop)

cut <- 11
drugSensitivity(matrix.MERGED, cell.line, cut)
plot( density( matrix.MERGED[,cell.line], na.rm=TRUE), main = "Full Set", xlab = "DSS" )
hist(matrix.MERGED[,cell.line])

drugs.sensitive <- topSensitive(matrix.MERGED, cell.line, cut)
drugs.resistant <- topResistant(matrix.MERGED, cell.line)

drug.list <- drugs.sensitive$DrugName
x <- matrix.MERGED[drug.list,]

# 3. Upload corresponding data set wit clusters
load('RData/leukemiaClust.RData')


# 4. Push Both Sets for Enrichment
# That is to verify that most of sensitive drugs 
# from a set tend to appear in the same cluster
enrichment.table <- buildEnrichmentD(tree.DRUGS, drugs.sensitive, drugs.resistant)

# Add information (new col) to 'tree.DRUGS' about
# which drugs to highlight: sensitive or resistant
is.top <- tree.DRUGS[,"DrugName"] %in% drugs.sensitive$DrugName # sensit
is.bot <- tree.DRUGS[,"DrugName"] %in% drugs.resistant$DrugName # resist
tree.DRUGS[,"isTop"] <- 0 
tree.DRUGS[is.top,"isTop"] <- 1
#tree.DRUGS[is.bot,"isTop"] <- -1

dropJSON(tree.DRUGS, path='/home/comrade/Projects/d3.v3/circular.json')


tree.drugs.TARGET <- data.frame( Cluster=integer(), PubChem.CID=factor(), DrugName=character(), 
                                    Kegg=factor(), Who=character(), Level1=factor(), Level1=factor() )
for(cluster in unique(tree.DRUGS$Cluster)){
  drop <- tree.DRUGS[,"Cluster"] == cluster & !is.na(tree.DRUGS[,"PubChem.CID"])
  cids <- tree.DRUGS[drop, c("DrugName","PubChem.CID")]
  for(cid in cids$PubChem.CID){
    if(cid %in% annotations.MERGED$Pubchem_CID){
      drop <- drop <- which( annotations.MERGED[,"Pubchem_CID"] == cid)
      keggs <- annotations.MERGED[drop,"KEGG_id"]
      for(kegg in keggs){
        who <- target.class.AML[target.class.AML[,"KEGG_id"] == kegg,"Who_Name"]
        l1 <- target.class.AML[target.class.AML[,"KEGG_id"] == kegg,"level_1_Description"]
        l2 <- target.class.AML[target.class.AML[,"KEGG_id"] == kegg,"level_2_Description"]
        df <- data.frame( Cluster=as.integer(cluster), PubChem.CID=cid, 
                            DrugName=cids[cids[,"PubChem.CID"] == cid,"DrugName"], 
                              Kegg=kegg, Who=who, Level1=factor(l1), Level2=factor(l2) )
        
        colnames(df) <- c("Cluster","PubChem.CID","DrugName","Kegg","Who","Level1", "Level2")
        tree.drugs.TARGET <- rbind(tree.drugs.TARGET, df)
      }       
    }
  }
}
remove(drop,who,l1,l2,df,cid,cids,kegg,keggs,cluster)
p <- tree.drugs.TARGET[tree.drugs.TARGET[,"Cluster"] == 8, c(5,6,7)]
plotDrugClassesDistibution(p, category.name='Targets')  

tab <- table(tree.drugs.TARGET$Cluster,tree.drugs.TARGET$Level1)
tab <- data.frame(tab)
colnames(tab) <- c("Cluster", "SampleName", "isTop")
tab$Cluster <- as.integer(tab$Cluster)


#tab.sep <- tab[tab[,"Cluster"] == 4 | tab[,"Cluster"] == 2,]
tab.sep <- tab
drop <- which( tab.sep$isTop == 0)
tab.sep <- tab.sep[-drop,]
dropCirclePackingJSON(tab.sep, path='/home/comrade/Projects/d3.v3/circle_packing.json')

