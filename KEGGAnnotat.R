#
# This is to sort out initial (raw) data
# and generate comprehensive FIMM Drug Dictionary
# along with Annotations for those drugs;
# HTML Summary included.
#
#######################################################

setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('pipeline_sup.R')
require(Nozzle.R1)

REPORT <- newCustomReport( "FIMM Collection:", asEmph(" Summary.") )

### Import Data Sets and Annotations ###
# ==================================== #

read.csv(file="../datasets/KEGG_ATC.tsv", head=TRUE, sep="\t") -> annotations.KEGG
#read.csv(file="../datasets/drug_dictionary_FIMM.csv", head=TRUE, sep="\t") -> dict.FIMM
load('RData/FimmCollection.RData')
collection.FIMM -> dict.FIMM
remove(collection.FIMM)

# Load four drug classifications from KEGG
read.csv(file="../datasets/1_biological_classes_kegg.csv", head=TRUE, sep=",") -> bio.class.FIMM
read.csv(file="../datasets/2_usp_classes_kegg.csv", head=TRUE, sep=",") -> usp.class.FIMM
read.csv(file="../datasets/3_target_classes_kegg.csv", head=TRUE, sep=",") -> target.class.FIMM
read.csv(file="../datasets/4_antineoplastics_classes_kegg.csv", head=TRUE, sep=",") -> neoplast.class.FIMM


### Produce FIMM Annotations  ###
# ============================= #
fimmcollection.section.REPORT <- newSection( "Fimm Collection Annotations" )

fimm.drug.count <- length( unique(dict.FIMM$FIMM.Batch.ID) )

# Drop drugs without PubChem.CID in FIMM Dctionary
drop <- which( is.na(dict.FIMM$PubChem.CID) )
dict.FIMM <- dict.FIMM[-drop,]

# Add if a drug has KEGG.ID
add <- dict.FIMM[,"PubChem.CID"] %in% unique(annotations.KEGG[,"Pubchem_CID"])
dict.FIMM$hasKegg <- add

# Count unique PUBCHEM IDs
fimm.drug.CID <- length( unique(dict.FIMM$PubChem.CID) )
summary.FIMM <- data.frame(FIMM.Batch=fimm.drug.count, PubChem.CIDs=fimm.drug.CID)
remove(fimm.drug.CID, fimm.drug.count)  # clean up

# Build annotation table for FIMM Collection:
# annotation.FIMM supposed to link PubChem.CID with KEGG_ID
drop <- annotations.KEGG[,"Pubchem_CID"] %in% unique(dict.FIMM[,"PubChem.CID"]) 
annotations.FIMM <- annotations.KEGG[drop,-c(18,19,20)]
remove(annotations.KEGG)

annotations.FIMM <- unique(annotations.FIMM)  # remove duplicated rows

# Separate away ATC classes from annotation table
# "cols.to.extract" corresponds to the following colums:
# 'Who_Name' and Levels: 1 to 5 along with descriptions
cols.to.extract <- c(7,8,9,10,11,12,13,14,15,16,17)  
atc.classes.FIMM <- data.frame(annotations.FIMM[,c( cols.to.extract ,3 )])  # 3 is for KEGG_id 
annotations.FIMM <- annotations.FIMM[,-cols.to.extract]
remove(cols.to.extract)

summary.FIMM$Annotations <- length(annotations.FIMM$Pubchem_CID)
summary.FIMM$annotations.KEGG <- length( unique(annotations.FIMM[,"Pubchem_CID"]) )

# (1) Separate KEGG classes so that *.class.FIMM
# contain drugs (KEGG.IDs) only from Fimm collecrion
# (2) Each *.class.FIMM table has to be of the same length 
# as annotations.FIMM. It means that drugs with missing annotations
# are assigned to 'undefined' class.
# BIO
drop <- bio.class.FIMM$KEGG_id %in% annotations.FIMM$KEGG_id
bio.class.FIMM <- bio.class.FIMM[drop,]  # (1)
drop <- !(annotations.FIMM$KEGG_id %in% bio.class.FIMM$KEGG_id)
ids <- unique(annotations.FIMM[drop, "KEGG_id"])  # (2)
bio.class.FIMM <- fulfillClassFrame(bio.class.FIMM, ids)
# TARGET
drop <- target.class.FIMM$KEGG_id %in% annotations.FIMM$KEGG_id
target.class.FIMM <- target.class.FIMM[drop,]
drop <- !(annotations.FIMM$KEGG_id %in% target.class.FIMM$KEGG_id)
ids <- unique(annotations.FIMM[drop, "KEGG_id"])  # (2)
target.class.FIMM <- fulfillClassFrame(target.class.FIMM, ids)
# USP
drop <- usp.class.FIMM$KEGG_id %in% annotations.FIMM$KEGG_id
usp.class.FIMM <- usp.class.FIMM[drop,]
drop <- !(annotations.FIMM$KEGG_id %in% usp.class.FIMM$KEGG_id)
ids <- unique(annotations.FIMM[drop, "KEGG_id"])  # (2)
usp.class.FIMM <- fulfillClassFrame(usp.class.FIMM, ids)
# NEOPLAST
drop <- neoplast.class.FIMM$KEGG_id %in% annotations.FIMM$KEGG_id
neoplast.class.FIMM <- neoplast.class.FIMM[drop,]
drop <- !(annotations.FIMM$KEGG_id %in% neoplast.class.FIMM$KEGG_id)
ids <- unique(annotations.FIMM[drop, "KEGG_id"])  # (2)
neoplast.class.FIMM <- fulfillClassFrame(neoplast.class.FIMM, ids)

# clean up factor levels for classes frames
neoplast.class.FIMM <- cleanUpFactors(neoplast.class.FIMM)
target.class.FIMM <- cleanUpFactors(target.class.FIMM)
atc.classes.FIMM <- cleanUpFactors(atc.classes.FIMM)
bio.class.FIMM <- cleanUpFactors(bio.class.FIMM)
usp.class.FIMM <- cleanUpFactors(usp.class.FIMM)


plotDrugClassesDistibution(target.class.FIMM, category.name='Targets')
plotDrugClassesDistibution(neoplast.class.FIMM, category.name='Antineoplastic: ')
plotDrugClassesDistibution(atc.classes.FIMM, category.name='ATC: ')
plotDrugClassesDistibution(bio.class.FIMM, category.name='Biological')
plotDrugClassesDistibution(usp.class.FIMM, category.name='USP')

fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( summary.FIMM, "FIMM Collection Summary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( dict.FIMM[1:7,], "FIMM Dictionary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( annotations.FIMM[1:3,c(1,3,8,9,10,11)], "FIMM Annotations" ) )
remove(drop)

# Save R Objects to a file
save(annotations.FIMM, dict.FIMM, summary.FIMM, bio.class.FIMM, usp.class.FIMM, target.class.FIMM, neoplast.class.FIMM, file = "RData/FimmDrugAnnotations.RData")

REPORT <- addTo( REPORT, fimmcollection.section.REPORT );
writeReport( REPORT, filename="../reports/FIMMCollectionKEGG/REPORT" )