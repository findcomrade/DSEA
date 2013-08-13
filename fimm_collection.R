#
# This is to sort out initial (raw) data
# and generate comprehensive FIMM Drug Dictionary
# along with Annotations for those drugs;
# HTML Summary included.
#
#######################################################

setwd("/home/comrade/Ubuntu One/DSEA/r-code")
# load("../datasets/merged_chembl.RData")

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

fimm.drug.count <- length( unique(fimm.DICT$FIMM.batch.ID) )
# Drop drugs without pubchem_CIDs in FIMM Dctionary
drop <- which( is.na(fimm.DICT$pubchem_CIDs) )
fimm.DICT <- fimm.DICT[-drop,]

# Count unique PUBCHEM IDs
fimm.drug.CID <- length( unique(fimm.DICT$pubchem_CIDs) )
fimm.SUMMARY <- data.frame(FIMM.Batch=fimm.drug.count, PubChem.CIDs=fimm.drug.CID)
remove(fimm.drug.CID, fimm.drug.count)  # clean up

# Build annotation table for FIMM Collection
drop <- atc.ANNO[,"Pubchem_CID"] %in% unique(fimm.DICT[,"pubchem_CIDs"]) 
fimm.ANNO <- atc.ANNO[drop,-c(19,20)]

# Drop Rows (compounds) without Annotations
drop <- which( is.na(fimm.ANNO$Level1) )
missing.ANNO <- fimm.ANNO[drop,]
fimm.ANNO <- fimm.ANNO[-drop,]
fimm.ANNO <- unique(fimm.ANNO)  # remove duplicated rows
fimm.SUMMARY$Annotations <- length(fimm.ANNO$Pubchem_CID)
fimm.SUMMARY$PubChem.Annotated <- length( unique(fimm.ANNO[,"Pubchem_CID"]) )
# fimm.SUMMARY$Missing.Annotat <- length( unique(missing.ANNO[,"Pubchem_CID"]) )

fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.SUMMARY, "FIMM Collection Summary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.DICT[1:7,], "FIMM Dictionary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.ANNO[1:3,c(1,8,9,10,11)], "FIMM Annotations" ) )
remove(atc.ANNO, drop)

# Save R Objects to a file
save(fimm.ANNO, missing.ANNO, fimm.DICT, fimm.SUMMARY, file = "FimmATCDrugAnnotations.RData")

REPORT <- addTo( REPORT, fimmcollection.section.REPORT );
writeReport( REPORT, filename="../reports/FIMMCollection/REPORT" )