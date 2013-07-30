#
# MAIN FILE which implements DSEA 
# (Drug Set Enrichment Analysis) pipeline;
# 
# Author: Dmitrii Bychkov, FIMM 2013
# (dmitrii.bychkov@helsinki.fi)
#######################################################

setwd("/home/comrade/Ubuntu One/DSEA/r-code")
dir.create( "../reports", showWarnings=FALSE );

library(vcd)
library(grid)
library(lattice)
library(ggplot2)
require(Nozzle.R1)

source("DSEA_SUP.r")

# LOAD DATA
ds.name <- "sanger"
FRAME <- loadDataSet(ds.name)
cell.lines <- getCellLines(ds.name)
pd <- pData(sangerSet)

# VISUALIZE (SENSIVITY) for ONE of the CELL LINES
plotDrugSensitivity("T47D")  # TODO: implement choice of the cell.line

### DEFINITIONS:###
target.cell.line <- "ES5"
target.top.count <- length(getTopDrugs(target.cell.line))            # top drugs (choose according to the response measurement) 
threshold <- 0.05            # higher p-value
###################

# CALCULATE ENRICHMENT
# target.CL - the one we paste to get enrichment
# reference.CL - all those that are already available in a data set 
ENRICHMENT <- buildEnrichmentTable(target.cell.line)

# TARGET CELL LINE OVERVIEW
TARGET.CL.OVERVIEW <- pd[pd["Cell.Line"] == as.character(target.cell.line),]
colnames(TARGET.CL.OVERVIEW)[2] <- "top.count"
TARGET.CL.OVERVIEW$top.count <- target.top.count

# DIAGNOSIS ENRICHMENT
#plot(plotDiagnosisFreq()) # for initial data set
#plot( plotEnrichedDiagnos() )


# summary table for enrichment
SUMMARY <- ENRICHMENT[ENRICHMENT$significant,]
#SUMMARY$Percent <- paste( format(round(ENRICHMENT[ENRICHMENT$SCORE > threshold,"SCORE"]/target.top.count *100, 2), nsmall=2), " %", sep="" )


# REPORT: Create Report Elements №№№№№№№№№№№№№№№№№№№
REPORT <- newCustomReport( "Drug Set Enrichment Analysis (DSEA): <br>", asEmph(" A Computational Approach to Identify <br> Functional Drug Sets from High-Throughput Drug Testing.") )

# --- Introduction ---
introduction.SECTION <- newSection( "Introduction" )
introduction.SECTION <- addTo( introduction.SECTION, newParagraph(asStrong("BASICS:")) )
introduction.SECTION <- addTo( introduction.SECTION, newList(
                                                     newParagraph("High-Throughput Drug Testing allows to reveal the impact of a particular drug (set of multiple drugs) on an individual;"),    
                                                     newParagraph("While methodology of the biological experiment quite similar to <i>RNAi or Small Molecules Screening</i>..."),
                                                     newParagraph("...its computational side brings new challenges when identifying new drug sets.")    
                                                    ) 
                              )

introduction.SECTION <- addTo( introduction.SECTION, newParagraph(asStrong("HYPOTHESIS: "), " The most sensitive drugs in any given individual cancer sample tend to show similar response in a cohort of other similar samples.") )
figure.file <- paste("sens_example",".png",sep="")
png( paste( "../reports/", figure.file, sep="" ), width=750, height=500 )
plotDrugSensitivity("T47D") 
dev.off()
introduction.SECTION <- addTo( introduction.SECTION, newFigure(figure.file, "T47D Cell Line (Sanger Set: http://www.cancerrxgene.org/)") )
introduction.SECTION <- addTo( introduction.SECTION, newParagraph(asStrong("RELATION TO GENE SET ENRICHMENT ANALYSIS: "), " GSEA interprets gene expression data by focusing on gene sets, that is, groups of genes that share common biological function, chromosomal location, or regulation." ) )
introduction.SECTION <- addTo( introduction.SECTION, newParagraph(asStrong("Drug Set Enrichment Analysis intends to facilitate the assessment of drug responses across multiple cell lines and also allow for drug response prediction.")) ) 

# --- Data Set ---
dataset.SECTION <- newSection( "Data Set Overview" )
dataset.SECTION <- addTo( dataset.SECTION, newParagraph(asStrong("TOY DATA SET:")) )
dataset.SECTION <- addTo( dataset.SECTION, newParagraph("The Genomics of Drug Sensitivity in Cancer project - is an academic research program to identify molecular features of cancers that predict response to anti-cancer drugs.") )
dataset.SECTION <- addTo( dataset.SECTION, newParagraph(asSummary( newResult( asParameter( "Drugs (IC_50):  " ), " ", asValue( "139" ), isSignificant=FALSE ) )) )
dataset.SECTION <- addTo( dataset.SECTION, newParagraph(asSummary( newResult( asParameter( "Cell Lines:  " ), " ", asValue( "714" ), isSignificant=FALSE ) )) )
dataset.SECTION <- addTo( dataset.SECTION, newParagraph("see: http://www.cancerrxgene.org") )

dataset.SECTION <- addTo( dataset.SECTION, newParagraph(asStrong("CONTENT:")) )
dataset.overview <- cbind( gsub("1","",summary(pData(sangerSet))[,1]), summary(pData(sangerSet))[,3:4] )
colnames(dataset.overview)[1] <- "Cell.Line"
dataset.overview <- dataset.overview[-length(dataset.overview[,1]),]  # trim the lasr row
dataset.SECTION <- addTo( dataset.SECTION, newTable( dataset.overview, "Target Cell Line Overview" ) )

figure.file <- paste("sens_example",".png",sep="")
png( paste( "../reports/", figure.file, sep="" ), width=750, height=500 )
plotDrugSensitivity("T47D") 
dev.off()
dataset.SECTION <- addTo( dataset.SECTION, newFigure( figure.file, "T47D Cell Line" ) )

dataset.distr.SECTION <- newSection( "Tissue and Diagnosis Distribution" )
figure.file <- paste("tissue_freq",".png",sep="")
ggsave(plotTissueFreq(), filename=paste( "../reports/", figure.file, sep="" ), width=300, height=175, limitsize=FALSE, units="mm")
dataset.distr.SECTION <- addTo( dataset.distr.SECTION, newFigure( figure.file, "Tissue Frequency" ) )

figure.file <- paste("diagnosis_freq",".png",sep="")
ggsave(plotDiagnosisFreq(), filename=paste( "../reports/", figure.file, sep="" ), width=300, height=175, limitsize=FALSE, units="mm")
dataset.distr.SECTION <- addTo( dataset.distr.SECTION, newFigure( figure.file, "Diagnosis Frequency" ) )

dataset.SECTION <- addTo(dataset.SECTION, dataset.distr.SECTION)

# -- Methodology ---
methodology.SECTION <- newSection( "Methodology" )
methodology.SECTION <- addTo( methodology.SECTION, newParagraph(asStrong("PRELIMINARY STAGE: ")) )
methodology.SECTION <- addTo( methodology.SECTION, newList(
                                          newParagraph("Collect DSRT data from different sources (Need to bind drugs to unique IDs);"),    
                                          newParagraph("Cluster those drugs to identify functional drug sensitivity profiles;")
                                                            ) 
                        )

methodology.SECTION <- addTo( methodology.SECTION, newParagraph(asStrong("INPUT: "), " A new cancer case screen.") ) 

methodology.SECTION <- addTo( methodology.SECTION, newParagraph(asStrong("ENRICHMENT STAGE: ")) )
methodology.SECTION <- addTo( methodology.SECTION, newList( isNumbered=TRUE,
                                          newParagraph("Order drugs according to sensitivity;"),    
                                          newParagraph("Compute overlap;"),
                                          newParagraph("Apply Fisher's exact test to obtain enrichment p-value;"),
                                          newParagraph("Adjust for multiple testing;"),
                                          newParagraph("Choose cell lines with p-value < 0.05;"),
                                          newParagraph("Identify a set of prospectively sensitive drugs for a given sample (input);"),
                                          newParagraph("Access the enrichment of top-drugs from the functional drug classes.")
                                                            ) 
                        )


line1 <- target.cell.line
which <- 1
line2 <- as.character(ENRICHMENT[which,1])

# TARGET CELL LINE OVERVIEW
method.table1 <- rbind( pd[pd["Cell.Line"] == as.character(line1),], pd[pd["Cell.Line"] == as.character(line2),] )
colnames(method.table1)[2] <- "top.count"
method.table1$top.count <- c( length(getTopDrugs(line1)), length(getTopDrugs(line2)) )

methodology.SECTION <- addTo( methodology.SECTION, newTable(method.table1, paste(line1, " VS ", line2, sep="")) )


figure.file <- paste("method_example",".png",sep="")
figure.plot <- genOverlapPlot( target.cell.line, as.character(ENRICHMENT[which,1]), target.top.count, as.numeric(ENRICHMENT[which,2]) )
ggsave(figure.plot, filename=paste( "../reports/", figure.file, sep="" ), width=300, height=175, limitsize=FALSE, units="mm")
methodology.SECTION <- addTo( methodology.SECTION, newFigure( figure.file, "Overlapping Drugs Among Top Sensitive" ) )
                              
methodology.SECTION <- addTo( methodology.SECTION, newTable(ENRICHMENT[which,], paste(line1, " Enrichment ", sep="")) ) 
# (**) chose top according to p-value; also relevant to (*)
# for (cl in ENRICHMENT[ENRICHMENT$SCORE > threshold,1]){
#   figure.file <- paste(cl,".png",sep="")
#   
#   ggp <- genOverlapPlot(target.cell.line, cl, target.top.count)
#   ggsave(ggp, filename=paste( "reports/", figure.file, sep="" ), width=300, height=175, limitsize=FALSE, units="mm")
#   
#   section.VISUAL <- addTo( section.VISUAL, newFigure( figure.file, cl ) )
# }

# --- Enrichment ---
enrichment.SECTION <- newSection( "Results" )
enrichment.SECTION <- addTo( enrichment.SECTION, newParagraph( paste("Target Cell Line :: ",target.cell.line,sep="") ) )
enrichment.SECTION <- addTo( enrichment.SECTION, newTable( TARGET.CL.OVERVIEW, "Target Cell Line Overview" ) )
enrichment.SECTION <- addTo( enrichment.SECTION, newTable( SUMMARY, "Top Scored Cell Lines" ) )

figure.file <- paste("tissue_enrich",".png",sep="")
ggsave(plotEnrishedTissue(), filename=paste( "../reports/", figure.file, sep="" ), width=300, height=175, limitsize=FALSE, units="mm")
enrichment.SECTION <- addTo( enrichment.SECTION, newFigure( figure.file, "Enriched Tissues" ) )

figure.file <- paste("diagnosis_enrich",".png",sep="")
ggsave(plotEnrichedDiagnos(), filename=paste( "../reports/", figure.file, sep="" ), width=300, height=175, limitsize=FALSE, units="mm")
enrichment.SECTION <- addTo( enrichment.SECTION, newFigure( figure.file, "Enriched Diagnosis" ) )

# --- Acknowledgements ---
acknowl.SECTION <- newSection( "Acknowledgements" )
acknowl.SECTION <- addTo( acknowl.SECTION, newParagraph(asStrong("Supervisors:")) )
acknowl.SECTION <- addTo( acknowl.SECTION, newList( newParagraph( "MSc. John Patrick Mpindi" ), newParagraph( "Dr. Paivi Ostling" ), newParagraph( "Prof. Olli Kallioniemi" ) ) )

# REPORT: Assemble Report Structure Bottom-Up
REPORT <- addTo( REPORT, introduction.SECTION, dataset.SECTION, methodology.SECTION, enrichment.SECTION, acknowl.SECTION  );

# render report to file
writeReport( REPORT, filename="../reports/my_report" )
