setwd("/home/comrade/Ubuntu One/DSEA/r-code")
library(ggplot2)
library(grid)
library(XLConnect)
library(reshape)




DSS_table_wb <- loadWorkbook ("/home/comrade/Desktop/Merged_drug_screening_data.xlsx")
DSS_tbl <- readWorksheet(DSS_table_wb, sheet=1, header = TRUE)
#Filter out drugs screened over a few cell lines
filt <- apply(as.matrix(DSS_tbl[,-c(1:5)]), 1, function(x) sum(is.na(x)) < 0.50 * length(x))
DSS_tbl <- DSS_tbl[filt,]
DSS_tbl <- DSS_tbl[as.character(DSS_tbl$Name.Drug) %in% drug.list,]
drugs <- as.character(DSS_tbl$Name.Drug)
xpr <- as.matrix(DSS_tbl[,-c(1:5)])
rownames(xpr) <- drugs

xpr <- melt(xpr)
#fct <- cut(xpr$value, breaks = c(0,0.5,1 , max(xpr$value,na.rm=TRUE) * 2))
fct <- cut(xpr$value, breaks = c(min(xpr$value,na.rm=TRUE) * 2,5,15, max(xpr$value,na.rm=TRUE) * 2))
levels(fct) <- c("Resistant (< 5)","Intermediate (5 - 15)","Sensitive (>= 15)")
mat <-  table(xpr$X1,fct)
df <-  prop.table(mat,1)
df <- df[,3:1]
df <- df[order(df[,1],decreasing=FALSE),]
#Order of drug names. You can replace this order with the drug list you are querrying.
levels_order <- rownames(df)
df.m <- melt(df)
df.m <- rename(df.m, c(Var.1 = "Compound", fct = "Groups"))
df.m$Compound <- factor(df.m$Compound,levels=levels_order)
df.m$Groups <- factor(df.m$Groups,levels=c("Sensitive (>= 15)","Intermediate (5 - 15)","Resistant (< 5)"))

#Set the plot size
X11(width=10.67, height=10.69)
main_title <- paste("Sensitivity profile to",length(drugs),"compounds by DSS scores",sep=" ")
base_size=15
a <- ggplot(df.m, aes(x = Compound, y=value*100, fill = Groups)) + labs(x = "Compound", y = "Proportion out of 182 cancer samples (in percentage)") + opts(title = main_title,legend.key.size = unit(2, "lines"),
                                                                                                                                                           legend.text = theme_text(size = base_size * 0.7),legend.title = theme_text(size = base_size * 0.8, face = "bold", hjust = 0), axis.text.x = theme_text(angle = 90,hjust = 1, size = base_size),axis.title.x =  theme_text(size = base_size), axis.title.y = theme_text(size = base_size*0.9, angle = 90,vjust = 0.3),
                                                                                                                                                           panel.grid.major = theme_line(colour = "grey90"), axis.text.y = theme_text(hjust = 1,family= "mono",face= "bold",colour = "black", size = base_size*0.8),plot.title = theme_text(size = base_size),
                                                                                                                                                           panel.grid.minor = theme_blank(), panel.background= theme_blank(),plot.background = theme_rect(fill = "white", colour = "white" ))

#b <- a + geom_bar(stat = "identity", position = "stack")+scale_fill_discrete("Effective/not effective") + coord_flip()
b <- a + geom_bar(stat = "identity", position = "stack")+scale_fill_discrete("Response (DSS range)") + coord_flip()
b
ggsave(filename="/home/comrade/Desktop/SR_sensitive.png",dpi=300)

source('testing-related.r')

read.csv(file="../datasets/updated_22_aug/AML_updated.csv", head=TRUE, sep="\t") -> aml.DATA
read.csv(file="../datasets/drug_dictionary_FIMM.csv", head=TRUE, sep="\t") -> fimm.DICT
read.csv(file="../datasets/updated_22_aug/fimm_collection.csv", head=TRUE, sep="\t") -> fimm.UPDATE


fimm.UPDATE <- unique(fimm.UPDATE)


aml.DATA[,"ID.Drug"] <- strtrim(aml.DATA[,"ID.Drug"],10)  # trim batch.id
fimm.DICT$ID.Drug <- strtrim(fimm.DICT[,"FIMM.batch.ID"],10)

collection.FIMM <- cbind(fimm.DICT[,5], fimm.DICT[,1:4])
colnames(collection.FIMM) <- c("FIMM.ID", "FIMM.Batch.ID", "Drug.Name", "ChEMBL.ID", "PubChem.CID")

collection.FIMM <- unique(collection.FIMM)

sum( unique(fimm.UPDATE[,"FIMM.ID"]) %in% unique(collection.FIMM[,"FIMM.ID"]) )  # 256
sum(collection.FIMM[,"FIMM.ID"] %in% fimm.UPDATE[,"FIMM.ID"] )  # 256
sum( unique(fimm.UPDATE[,"FIMM.BATCH.ID"]) %in% unique(collection.FIMM[,"FIMM.Batch.ID"]) )  # 201
sum( unique(fimm.UPDATE[,"Drug.Name"]) %in% unique(collection.FIMM[,"Alias"]) )  # 242

# bind Targets to collection
collection.FIMM$Target <- "Undefined"
wh <- which( collection.FIMM[,"FIMM.ID"] %in% unique(fimm.UPDATE[,"FIMM.ID"]) )
for(row in wh){
  id <- as.character( collection.FIMM[row,"FIMM.ID"] )
  collection.FIMM[row, "Target"] <- as.character( unique( fimm.UPDATE[ fimm.UPDATE[,"FIMM.ID"] == id, "Mechanism.Targets"] ))
}

save(collection.FIMM, file = "RData/FimmCollection.RData")

sum( aml.DATA[,"ID.Drug"] %in% fimm.UPDATE[,"FIMM.ID"] ) 
sum( aml.DATA[,"ID.Drug"] %in% fimm.DICT[,"FIMM.batch.ID"] ) 

