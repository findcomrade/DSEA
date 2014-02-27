library(ggplot2)
library(XLConnect)
library(reshape)
library(grid)
library(Hmisc)

setwd("/home/comrade/Ubuntu One/DSEA/r-code")


# <----------  Inputs:        ----------> 
source.xlsx             <- "/home/comrade/Desktop/Merged_drug_screening_data.xlsx"
current_sample          <- "X564_23012012_9999_BM"
drug.feature            <- "Name.Drug"
top.drugs               <- 79                      # [10:99]
top.samples             <- 25                      # [1:39]
features.head           <- 2
# <----------  end of Inputs  ----------> 

DSS_table_wb  <- loadWorkbook (source.xlsx)
DSS_tbl       <- readWorksheet(DSS_table_wb, sheet=1, header = TRUE)
remove(DSS_table_wb)

top.samples   <- top.samples + 1
samples.count <- dim(DSS_tbl)[2] - features.head

#Filter out drugs screened over a few cell lines
#filt <- apply(as.matrix(DSS_tbl[,-c(1:5)]), 1, function(x) sum(is.na(x)) < 0.50 * length(x))

# Order drugs according to sensitivity in the sample
filt      <- order(DSS_tbl[,current_sample], decreasing=TRUE)
DSS_tbl   <- DSS_tbl[filt,]

# Trim top sensitive drugs
DSS_tbl   <- DSS_tbl[1:top.drugs,]

drugs <- as.character(DSS_tbl[,drug.feature])
xpr <- as.matrix(DSS_tbl[,-c(1:features.head)])
rownames(xpr) <- drugs

#tmp.enrich <- as.character( correlation.table[ ,"Cell.Line"] )
#tmp.enrich <- as.character( enrichment.table[enrichment.table[,"significant"] ,"Cell.Line"] )

#col.mat <- which( colnames(matrix.NEW) %in% tmp.enrich)
#colnames(matrix.NEW)[col.mat] <- strtrim(colnames(matrix.NEW)[col.mat], 18)
#tmp.enrich <- strtrim(tmp.enrich, 18)
#current_sample <- strtrim(current_sample, 18)
 
#indx <- which( rownames(matrix.NEW) %in% drugs.sensitive$DrugName )
#xpr <- matrix.NEW[indx,]
#drugs <- as.character(data.MERGED$Name.Drug)  

# Create Groups for Percentage Plot
# ==================================
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

# Now create dataframe of the single sample,
# aka Top Correlating Samples
# =========================================
# xpr_patients <- as.matrix(DSS_tbl[,-c(1:5)])
#xpr_patients <- xpr_patients[,idx_correlatedsamples]
#tmp.enrich <- as.character( correlation.table[ ,"Cell.Line"] )

#tmp.enrich <- as.character( enrichment.table[enrichment.table[,"significant"] ,"Cell.Line"] )
#tmp.enrich <- c( current_sample, tmp.enrich[1:20] )
#idx_correlatedsamples <- which( colnames(matrix.NEW) %in% tmp.enrich)

#xpr_patients <- matrix.NEW[,idx_correlatedsamples]

# Convert data table to numetic matrix
tmp.matrix <- as.matrix(DSS_tbl[,-c(1:features.head)])
rownames(tmp.matrix) <- drugs

# Compute correlations
corr.data <- rcorr(tmp.matrix)
rs <- corr.data$r[current_sample,]
ps <- corr.data$P[current_sample,]
corr.table <- data.frame(Corr=rs, pVal=ps)

# order according to correlation value
corr.table <- corr.table[ order(corr.table$Corr, decreasing=TRUE), ]
corr.table$Cell.Line <- as.factor(rownames(corr.table))

tmp.enrich <- corr.table[1:top.samples,"Cell.Line"]
idx_correlatedsamples <- which( colnames(tmp.matrix ) %in% tmp.enrich)
xpr_patients <- as.matrix(tmp.matrix[,idx_correlatedsamples])


#rownames(xpr_patients) <- drugs
xpr_patients <- melt(xpr_patients)

xpr_patients$X1 <- factor(as.character( xpr_patients$X1 ),levels=levels_order)
xpr_patients$X2 <- factor(as.character( xpr_patients$X2 ),levels=tmp.enrich)

xpr_patients <- xpr_patients[!is.na(xpr_patients[,"X1"]),]
xpr_patients <- data.frame(xpr_patients, data_type=factor( rep("Top correlating samples", nrow(xpr_patients)) ))

# Sensitivity Probability
# ================================
smp_xpr <- xpr[as.character(xpr$X2) %in% current_sample,]
colnames(smp_xpr) <- colnames(df.m)
smp_xpr$Groups <- as.character(smp_xpr$Groups)
smp_xpr$data_type <- rep(paste(current_sample,"Orig.DSS",sep=" "),nrow(smp_xpr))
smp_xpr$value <- ifelse(smp_xpr$value>0,smp_xpr$value*-1,smp_xpr$value)
df.m$data_type <- rep("Sensitivity probability",nrow(df.m))

# Merge the data frames
df.m$value <- df.m$value*100
df.m$data_type <- factor(df.m$data_type)
smp_xpr$data_type <- factor(smp_xpr$data_type)

#Set the plot size
#X11(width=18, height=10.6)
png(file = "Sanger_Sensitivity_Histogram_300dpi.png", width = 1750, height = 1300,units = "px",res=90)
main_title <- paste("Sensitivity profile to",length(drugs),"compounds by DSS scores",sep=" ")
base_size=15
a <- ggplot()+ labs( x = "Compound" , y = paste("DSS score ::: Percentage of \n sensitivity out of", samples.count, "samples", sep=" ") ) + theme(legend.key.size = unit(2, "lines"),legend.position=c(1,0.3),legend.justification=c(-0.1,0),plot.margin = unit(c(-0.7,6,0.5,0.5),"cm"),
legend.text = element_text(size = base_size * 0.7),legend.title = element_text(size = base_size * 0.8, face = "bold", hjust = 0), axis.text.x = element_text(angle = 51,hjust = 1, vjust = 1,family= "mono",face= "bold",colour = "black",size = base_size*0.85),axis.title.x =  element_text(size = base_size), axis.title.y = element_text(size = base_size*0.9, angle = 90,vjust = 0.5),
panel.grid.major = element_line(colour = "grey90"), axis.text.y = element_text(hjust = 1,family= "mono",face= "bold",colour = "black", size = base_size),plot.title = element_text(size = base_size), panel.grid.minor = element_blank(), panel.background= element_blank(),plot.background = element_rect(fill = "white", colour = "white" )) 
b <- a + geom_bar(data = df.m, aes(x = Compound, y=value, fill = Groups),stat = "identity", position = "stack") 
b <- b + geom_bar(data = smp_xpr, aes(x = Compound, y=value, fill = Groups),stat = "identity",position = "stack") 
b <- b + scale_fill_discrete("Response (DSS range)") 
b <- b + facet_grid(data_type ~. ,scales="free_y")+ scale_y_continuous(breaks=c(-100,-75,-50,-25,0,25,50,75,100),labels=c(100,75,50,25,0,"25%","50%","75%","100%"))
#b <- b 
#b

library(gridExtra)
b2 <- ggplot(data = xpr_patients, aes(x = X1, y=X2, group=value))+ geom_tile(aes(fill = value)) + theme(legend.key.size = unit(2, "lines"),axis.text.x=element_blank(),legend.position=c(1,0.3),legend.justification=c(-0.2,0),
legend.text = element_text(size = base_size * 0.7),legend.title = element_text(size = base_size * 0.8, face = "bold", hjust = 0), axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.3,family= "mono",face= "bold",colour = "black",size = base_size),axis.title.x =  element_text(size = base_size), axis.title.y = element_text(size = base_size*0.9, angle = 90,vjust = 0.5),
panel.grid.major = element_line(colour = "grey90"), axis.text.y = element_text(hjust = 1,family= "mono",face= "bold",colour = "black", size = base_size),plot.title = element_text(size = base_size), panel.grid.minor = element_blank(), panel.background= element_blank(),plot.background = element_rect(fill = "white", colour = "white" ))
b2 <- b2 + labs(title = main_title)
b2 <- b2 + labs(x = NULL, y = "Correlated samples") 
b2 <- b2 + theme(plot.margin = unit(c(0.5,6,0.5,0.5),"cm"),axis.ticks.x=element_blank(),axis.text.x=element_blank())
b2 <- b2 + geom_text(aes(fill = xpr_patients$value,size=0.5, stat = "identity",label = round(xpr_patients$value, 0))) + scale_size(range=c(2.5))
b2 <- b2 + scale_fill_gradient("Response (DSS values)",low = "white", high = "red")
#b2
b2 <- b2 + facet_grid(data_type ~. ,scales="free_y") 

gp1<- ggplot_gtable(ggplot_build(b))
gp2<- ggplot_gtable(ggplot_build(b2))
maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
#Set each to the maximum width
gp1$widths[2:3] <- maxWidth
gp2$widths[2:3] <- maxWidth
p1 <- grid.arrange(gp2,gp1,nrow=2,widths=c(4,4))
p1
dev.off()


#ggsave(filename="./Sanger_Sensitivity_Histogram_300dpi.png",p1,dpi=300)

