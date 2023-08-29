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

if (length(args)!=4) {
  stop("Call: Rscript epitope_landscape.R <escape data file> <antibody mapping file> <fold resistance file> <output directory>", call.=FALSE)
} else  {
  inputfile_escape_data <- args[1]
  inputfile_antibodymapping <- args[2]
  inputfile_foldresistance <- args[3]
  outputdir <- args[4]
}

#Paper input:
# inputfile_escape_data <- "DATA/escape_data.csv"
# inputfile_antibodymapping <- "DATA/antibody_classes.csv"
# inputfile_foldresistance <- "DATA/Fold_Resistance_DMS_Sites_Epitopes.csv"
# outputdir <- "results_epitopelandscape"


## Outputfiles:
outputfile_heatmapplot <- "escape_fractions_mean_per_greaneyclass.pdf"
outputfile_heatmapplot_mergeE21 <- "escape_fractions_mean_per_greaneyclass_aggregateE1E2.pdf"
outputfile_heatmapplot_FR <- "foldresistance_DMS_sites_epitopes.pdf"
outputfile_heatmap_EF_FR <- "Heatmap_EF_FR.pdf"

col_escape_fraction <- "mut_escape"
dir.create(outputdir, showWarnings = FALSE)


##################### READ IN escape data from bloom lab
escape_data <- read.csv(inputfile_escape_data)
antibodymapping <- read.csv(inputfile_antibodymapping)

## filter the inputfiles: remove serum, keep antibodies
escape_data <- escape_data[-which(escape_data$condition_type == "serum"),]

## merge escape_data.csv with antibody classes from escape_data_site.csv which is in antibody_classes.csv
escape_data <- merge(escape_data, antibodymapping, by = "condition" )

## escape_fraction = mut_escape and set max value to 1:
x <- which(escape_data[col_escape_fraction]>=0.99)
if (length(x) > 0){ escape_data[x,col_escape_fraction] <- 0.99}
rm(x)
write.csv(escape_data, file = paste0(outputdir,"/escape_data_merged_greaneyclasses.csv"),row.names = FALSE, quote = FALSE)


###################### Define RBD sites, antibodies, and antibody classes:
RBDsites <- sort(unique(escape_data$site))
print(paste("Number of RBD sites:", length(RBDsites)))
groups <- sort(unique(escape_data$group))
print(paste("Antibody classes: ",paste(groups,collapse = ",")))
antibodies <- sort(unique(escape_data$condition))
print(paste("Number of antibodies: ", length(antibodies)))


antibodymapping2 <- escape_data[,c("condition","group")]
antibodymapping2 <- antibodymapping2[!is.na(antibodymapping2$group),]
antibodymapping2 <- antibodymapping2[!duplicated(antibodymapping2),]
print("Number of antibodies per group:")
table(antibodymapping2$group)
write.csv(antibodymapping2, file = paste0(outputdir,"/antibodymapping_greaneyclasses.csv"),row.names = FALSE, quote = FALSE)


# merge groups E1, E2.1 and E2.2 into one group
escape_data$group = ifelse(grepl("E2", escape_data$group, fixed = TRUE), "E12", escape_data$group)
escape_data$group = ifelse(escape_data$group == "E1", "E12", escape_data$group)

# aggregate mut_escape to site level
escape_data_site = aggregate(mut_escape ~ condition + site + group + IC50,
                             data = escape_data,
                             FUN = mean)

write.csv(escape_data_site, file = paste0(outputdir,"/dms_per_ab_per_site.csv"),row.names = FALSE, quote = FALSE)

####################### compare the AB classes by their mean escape fraction per site over all antibodies in this class.
EL <- as.data.frame(matrix(nrow = length(groups), ncol = length(RBDsites)))
colnames(EL) <- RBDsites; rownames(EL) <- groups

for (i in 1:length(groups)){
    g <- groups[i]
    #all sites in group i:
    g_escape_data <- escape_data[which(escape_data$group == g),]
    #subset by group g and aggregate over the mean mut escape at each site:
    a_mean <- aggregate(g_escape_data$mut_escape, list(g_escape_data$site), FUN=mean)
    EL[i,unlist(lapply(a_mean$Group.1,function(x) which(colnames(EL) == x)))] <- a_mean$x
}
rm(g);rm(g_escape_data);rm(a_mean)
#write.csv(EL, file=paste0(outputdir,"/escape_fractions_mean_per_greaneyclass.tsv"))
  
#remove sites with empty entries only:
j <- c()
for (i in 1:ncol(EL)){
    if (length(which(!is.na(EL[,i]))) == 0){
      j <- c(j, i)
    }
}
if (length(j) > 0) {
    EL2 <- EL[,-j]
}else{
    EL2 <- EL
}
EL2 <- EL2[,order(colnames(EL2))]
rm(j)
write.csv(EL2, file=paste0(outputdir,"/escape_fractions_mean_per_greaneyclass.tsv"))



######################### Merge E2 subclasses with E1

EL3 <- EL2
x <- unlist(lapply(c("E1","E2.1","E2.2"), function(x) which(rownames(EL3)==x)))
e12 <- apply(EL2[x,], 2 , function(x) mean(x, na.rm = T))
EL3 <- EL3[-x,]
EL3 <- rbind(EL3,e12 )
rownames(EL3)[nrow(EL3)] <- "E12"
EL3 <- EL3[order(rownames(EL3)),]
rm(x)



######################### Heatmap Plotting

pdf(paste0(outputdir,"/",outputfile_heatmapplot_mergeE21), width = 20)
heatmap.2(as.matrix(EL3), main = paste0("Escape fractions (mean) per site"), col=c("white",brewer.pal(8,"Blues")), symm = T,
          density.info = "none", trace="none",  
          Colv = FALSE, Rowv = FALSE, dendrogram = "none", scale = "none",
          margins = c(5, 5), sepwidth = c(0.05,0.05),
          sepcolor = "lightgrey",colsep = 1:ncol(EL2),rowsep = 1:nrow(EL2),
          na.color="grey",
          lhei=c(2, 10), cexRow = 1, cexCol = 0.5)
dev.off()

EL3lg <- log(EL3)
pdf(paste0(outputdir,"/",str_replace(outputfile_heatmapplot_mergeE21,".pdf","_log.pdf")), width = 20)
heatmap.2(as.matrix(EL3lg), main = paste0("Escape fractions (mean) per site, log scale"), col=c(brewer.pal(8,"Blues")), 
          symm = F, symkey = F, symbreaks = F,
          density.info = "none", trace="none",  Colv = FALSE, Rowv = FALSE, dendrogram = "none", scale = "none",
          margins = c(5, 5), sepwidth = c(0.05,0.05),sepcolor = "lightgrey",colsep = 1:ncol(EL2),rowsep = 1:nrow(EL2),
          na.color="grey",
          lhei=c(2, 10), cexRow = 1, cexCol = 0.5)
dev.off()
write.csv(EL3, file=paste0(outputdir,"/escape_fractions_mean_per_greaneyclass_aggregateE1E2.tsv"))




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

#Remove NTD:
# FRlg2_nontd <- FRlg2[-11,which(colnames(FRlg2)>318)]
# 
# pdf(paste0(outputdir,"/",outputfile_heatmapplot_FR), width = 20)
# heatmap.2(as.matrix(FRlg2_nontd), main = paste0("Fold Resistance per site (log 10 scale) DMS"), col=c("white",brewer.pal(8,"Blues")[2:8],"#053061", "black"), symm = T,
#           density.info = "none", trace="none",  Colv = FALSE, Rowv = FALSE, dendrogram = "none", scale = "none",
#           #margins = c(5, 5), sepwidth = c(0.05,0.05),sepcolor = "lightgrey",colsep = 1:ncol(FR),rowsep = 1:nrow(FR),
#           na.color="grey",
#           lhei=c(2, 10), cexRow = 1, cexCol = 0.5)
# dev.off()

