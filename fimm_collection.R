setwd("/home/comrade/Ubuntu One/DSEA/r-code")
#load("../datasets/merged_chembl.RData")

library(vcd)
library(grid)
library(gplots)
library(lattice)
library(ggplot2)
library(reshape2)
require(Nozzle.R1)

REPORT <- newCustomReport( "FIMM Collection:", asEmph(" Summary.") )
### Import Data Sets and Annotations ###
# ==================================== #
# merged_chembl_ATC_SCID -> atc.ANNO  # Chembl Drug Annotations
read.csv(file="../datasets/KEGG_ATC.tsv", head=TRUE, sep="\t") -> atc.ANNO
read.csv(file="../datasets/drug_dictionary_FIMM.csv", head=TRUE, sep=",") -> fimm.DICT
remove(merged_chembl_ATC_SCID)  # cleanup

### Produce FIMM Annotations  ###
# ============================= #
fimmcollection.section.REPORT <- newSection( "Fimm Collection Annotations" )

# Extract relevant information from 'ChEMBL db'
# to annotate drugs in FIMM data base.
fimm.drug.count <- length( unique(fimm.DICT$FIMM.batch.ID) )
drop <- which(is.na(fimm.DICT$pubchem_CIDs))
fimm.DICT <- fimm.DICT[-drop,]  # drop drugs without pubchem_CIDs in FIMM Dctionary
fimm.drug.CID <- length( unique(fimm.DICT$pubchem_CIDs) )
fimm.drug.CHE <- length( unique(fimm.DICT$CHEMBL_ID) )
fimm.SUMMARY <- data.frame(FIMM.Batch=fimm.drug.count, PubChem.CIDs=fimm.drug.CID, Chembl.IDs=fimm.drug.CHE)
remove(fimm.drug.CID, fimm.drug.count, fimm.drug.CHE)

# Build annotation table for FIMM Collection
drop <- atc.ANNO[,"Pubchem_CID"] %in% unique(fimm.DICT[,"pubchem_CIDs"]) 
fimm.ANNO <- atc.ANNO[drop,-c(19,20)]

# drop rows without annotations
drop <- which(is.na(fimm.ANNO$Level1))
missing.ANNO <- fimm.ANNO[drop,]
fimm.ANNO <- fimm.ANNO[-drop,]
fimm.ANNO <- unique(fimm.ANNO)  # remove duplicated rows
fimm.SUMMARY$PubChem.Annotated <- length(fimm.ANNO$Pubchem_CID)

fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.SUMMARY, "FIMM Collection Summary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.DICT[1:7,], "FIMM Dictionary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.ANNO[1:3,c(1,8,9,10,11)], "FIMM Annotations" ) )
remove(atc.ANNO, drop)


REPORT <- addTo( REPORT, fimmcollection.section.REPORT );
writeReport( REPORT, filename="../reports/FIMMCollection/REPORT" )