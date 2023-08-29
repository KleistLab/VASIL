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

if (length(args)!=3) {
  stop("Call: Rscript epitope_landscape.R <escape data file> <antibody mapping file> <output directory>", call.=FALSE)
} else  {
  inputfile_escape_data <- args[1]
  inputfile_antibodymapping <- args[2]
  outputdir <- args[3]
}


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
                             FUN = mean,
			     na.rm = TRUE)

write.csv(escape_data_site, file = paste0(outputdir,"/dms_per_ab_per_site.csv"),row.names = FALSE, quote = FALSE)



