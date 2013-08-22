#
# This is to sort out initial (raw) data
# and generate comprehensive FIMM Drug Dictionary
# along with Annotations for those drugs;
# HTML Summary included.
#
#######################################################

setwd("/home/comrade/Ubuntu One/DSEA/r-code")
require(Nozzle.R1)

REPORT <- newCustomReport( "FIMM Collection:", asEmph(" Summary.") )

### Import Data Sets and Annotations ###
# ==================================== #

read.csv(file="../datasets/KEGG_ATC.tsv", head=TRUE, sep="\t") -> kegg.ANNO
read.csv(file="../datasets/drug_dictionary_FIMM.csv", head=TRUE, sep=",") -> fimm.DICT

# Load four drug classifications from KEGG
read.csv(file="../datasets/1_biological_classes_kegg.csv", head=TRUE, sep=",") -> bio.KEGG
read.csv(file="../datasets/2_usp_classes_kegg.csv", head=TRUE, sep=",") -> usp.KEGG
read.csv(file="../datasets/3_target_classes_kegg.csv", head=TRUE, sep=",") -> target.KEGG
read.csv(file="../datasets/4_antineoplastics_classes_kegg.csv", head=TRUE, sep=",") -> neo.KEGG


### Produce FIMM Annotations  ###
# ============================= #
fimmcollection.section.REPORT <- newSection( "Fimm Collection Annotations" )

fimm.drug.count <- length( unique(fimm.DICT$FIMM.batch.ID) )
# Drop drugs without pubchem_CIDs in FIMM Dctionary
drop <- which( is.na(fimm.DICT$pubchem_CIDs) )
fimm.DICT <- fimm.DICT[-drop,]

# Add KEGG ID to Fimm dictionary
add <- fimm.DICT[,"pubchem_CIDs"] %in% unique(kegg.ANNO[,"Pubchem_CID"])
fimm.DICT$hasKegg <- add

# Count unique PUBCHEM IDs
fimm.drug.CID <- length( unique(fimm.DICT$pubchem_CIDs) )
fimm.SUMMARY <- data.frame(FIMM.Batch=fimm.drug.count, PubChem.CIDs=fimm.drug.CID)
remove(fimm.drug.CID, fimm.drug.count)  # clean up

# Build annotation table for FIMM Collection
drop <- kegg.ANNO[,"Pubchem_CID"] %in% unique(fimm.DICT[,"pubchem_CIDs"]) 
fimm.ANNO <- kegg.ANNO[drop,-c(18,19,20)]

fimm.ANNO <- unique(fimm.ANNO)  # remove duplicated rows
fimm.SUMMARY$Annotations <- length(fimm.ANNO$Pubchem_CID)
fimm.SUMMARY$PubChem.Annotated <- length( unique(fimm.ANNO[,"Pubchem_CID"]) )


fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.SUMMARY, "FIMM Collection Summary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.DICT[1:7,], "FIMM Dictionary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.ANNO[1:3,c(1,3,8,9,10,11)], "FIMM Annotations" ) )
remove(drop)

# Save R Objects to a file
save(fimm.ANNO, fimm.DICT, fimm.SUMMARY, file = "FimmKEGGDrugAnnotations.RData")

REPORT <- addTo( REPORT, fimmcollection.section.REPORT );
writeReport( REPORT, filename="../reports/FIMMCollectionKEGG/REPORT" )