#!/usr/bin/env Rscript

# Uses escape values from escape_data.csv to compute epitope landscape for all antibodies in a class, as well as across 
# all antibodies in a class. 
# For aggregating antibodies per class at each site, their mean values are being used. 
# Classes E1, E2.1 and E2.2 are being merged into E12


library(stringr)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(readr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Call: Rscript FR_sites_landscape.R <fold resistance file> <output dir>", call.=FALSE)
} else  {
  inputfile_foldresistance <- args[1]
  outputdir <- args[2]
}

## Outputfiles:
outputfile_heatmapplot_FR <- "foldresistance_DMS_sites_epitopes.pdf"

######################### FR Map:
FR <- read.csv(inputfile_foldresistance)
rownames(FR) <- FR$Epitope.Classes
FR <- FR[,-c(1,2)]
colnames(FR) <- str_replace(colnames(FR),"X","")

#Limit FR to 100 and log:
FR2 <- FR
FR2[FR2>100] <- 100
FRlg2 <- log10(FR2)

pdf(paste0(outputdir,"/",outputfile_heatmapplot_FR), width = 20)
heatmap.2(as.matrix(FRlg2), main = paste0("Fold Resistance per site (log 10 scale) DMS"), col=c("white",brewer.pal(8,"Blues")[2:8],"#053061", "black"), symm = T,
          density.info = "none", trace="none",  
          Colv = FALSE, Rowv = FALSE, dendrogram = "none", scale = "none",
          na.color="grey",
          lhei=c(2, 10), cexRow = 1, cexCol = 0.5)
dev.off()
