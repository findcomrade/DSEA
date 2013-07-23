setwd("/home/comrade/Ubuntu One/DSEA/r-code")
load("../datasets/merged_chembl.RData")
load("../datasets/sanger.RData")
load("../datasets/oregon.RData")

### Import Data Sets and Annotations ###
# ==================================== #
merged_chembl_ATC_SCID -> chembl.ANNO  # Chembl Drug Annotations

data.frame(exprs(sangerSet)) -> sanger.DATA
data.frame(exprs(xprOregonSet)) -> oregon.DATA

read.csv(file="../datasets/all_leukemia_cl_june_13_disha_astrid.csv", head=TRUE, sep=",") -> leukemia.DATA
read.csv(file="../datasets/original_dss_june_13_ovarian_akira_astrid.csv", head=TRUE, sep=",") -> ovarian.DATA
read.csv(file="../datasets/original_dss_breast_cancer_tamoxifen_project_sk_sh.csv", head=TRUE, sep=",") -> breast.DATA

read.csv(file="../datasets/fimm_chembl_dict.csv", head=TRUE, sep=",") -> fimm.DICT
read.csv(file="../datasets/sanger_chembl_dict.csv", head=TRUE, sep=",") -> sanger.DICT
read.csv(file="../datasets/oregon_chembl_dict.csv", head=TRUE, sep=",") -> oregon.DICT

pData(sangerSet) -> sanger.sample.ANNO
pData(xprOregonSet) -> oregon.sample.ANNO

# cleanup
remove(sangerSet, xprOregonSet, merged_chembl_ATC_SCID)


### Small Fixes for Initial Data Sets ###
# ===================================== #

# Snager Set:
sanger.DATA <- sanger.DATA[1:138,]  # Take only the first 138 drugs (those with IC_50 values)
rownames(sanger.DATA) <- gsub("_IC_50$", "", rownames(sanger.DATA))  # Trim "_IC_50$" from drug names 

# Oregon Set:
colnames(oregon.DATA) <- rownames(oregon.sample.ANNO)  # Trim "X" from sample name

# Ovarian Set:
# assign "NA" to empty FIMM_ID fields (there are 8 empty fields)
ovarian.DATA[-grep( "FIMM", as.character(ovarian.DATA[,1]), ignore.case = FALSE ), 1] <- NA




