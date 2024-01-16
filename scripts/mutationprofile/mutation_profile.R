#!/usr/bin/env Rscript

library(stringr)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(readr)
library(pheatmap)

## Generates mutation profile for set of variants from covsonar.
## Takes input data as obtained by covsonar (each lineage one tsv file) , 
## filters for mutations with at least 75% prevalence in all samples per lineage,
## greps for spike proteins, RBD sites. 
##
### INPUT: 
#### (In1) input data file as tsv as obtained by covsonar 
#### (In2) results directory
#### (In3) prefix for output files
#### (In4) threshold: mutations should be present in at least x% of all samples
#### (In5) date start: date to starting simulations
#### (In6) date end: date to end simulations
#### (In7) dms_per_ab_per_site: contains all positions that have dms data
### OUTPUT:
#### (Out1) <prefix>_mutations_spike_lists.csv : table with mutated sites for each lineage
#### (Out2) <prefix>_mutations_spike.csv : matrix with mutation profile for each lineage
#### (Out3) <prefix>_RBD_NTD_mutations.csv : table with mutational profile for each lineage.
#### (Out4) <prefix>_RBD_NTD_pseudogroups.csv: table with lineages and grouping based on identical mutational profile in RBD and NTD
#### (Out5) <prefix>_mutations_spikenumber_of_genomes_per_lineage.csv : table with number of genomes per lineage
#### (Out6) <prefix>_positiongroups_RBD_NTD_groups.pdf : Mutation Profile as plot
#### (Out7) <prefix>_positiongroups_RBD_NTD_groups_zoom.pdf : Mutation Profile of predefined lineages.

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop("Call: Rscript generate_mutation_profile <covsonar data file> <output directory> <prefix> <mutation threshold> <date start> <date end> <dms data>", call.=FALSE)
} else  {
  input_datafile_covsonar <- args[1]
  outputdir <- args[2]
  output_prefix_for_file <- args[3]
  threshold <- as.numeric(args[4])
  date_start <- args[5]
  date_end <- args[6]
  dms_per_ab_per_site_file <- args[7]
}


#Paper input:
dir.create(outputdir, showWarnings = FALSE)

## Outputfiles:
outputfile_mutationlist <- paste0(output_prefix_for_file,"_mutations_spike_lists.csv")
outputfile_mutationprofile <- paste0(output_prefix_for_file,"_mutations_spike.csv")
outputfile_mutationprofile_mutations <- paste0(output_prefix_for_file,"_RBD_NTD_mutations.csv")
outputfile_mutationprofile_plot <- paste0(output_prefix_for_file,"_positiongroups_RBD_NTD_groups.pdf")
outputfile_mutationprofile_groups <- paste0(output_prefix_for_file,"_RBD_NTD_pseudogroups.csv")




######## CONSTANT DEFINTIONS:
## spike protein positions (NTD and RBD):
#NTD_position_start <- 13
#NTD_position_end <- 317
#RBD_position_start <- 331
#RBD_position_end <- 531
dms_per_ab_per_site <- read.csv(dms_per_ab_per_site_file, sep = ",")
RBD_positions <- sort(unique(dms_per_ab_per_site$site))
RBD_position_start = min(RBD_positions)
RBD_position_end = max(RBD_positions)
NTD_positions <- c(14:20, 140:158, 245:264)
NTD_position_start = min(NTD_positions)
NTD_position_end = max(NTD_positions)


######################################### DATA PROCESSING 1: mutation profile

#### READ DATA and FILTER for mutations with at least 75% prevalence in all samples per lineage
D <- read.csv(input_datafile_covsonar, sep = "\t")
# Produce demo file
#v1 <- sort(sample(11:length(D$date)-10, 200, replace=F))
#t_keep <- c(1:10, v1)
#t_keep <- c(t_keep, (length(D$date)-9):length(D$date))
#print(length(t_keep))
#Dmock <- D[t_keep,]
#write.table(Dmock, file='covsonar_mock.tsv', quote=FALSE, sep='\t', col.names = NA)

### order dates
D <- D[order(D$date),]

## If actual data frame is smaller than start or end date:
#date_start <- max(date_start, min(D$date))
#date_end <- min(date_end, max(D$date))
if (!(date_start %in% D$date)){date_start = D$date[1]} ## indexes start at 1 in R not 0
if (!(date_end %in% D$date)){date_end = D$date[length(D$date)]}



## Check dates:
dates <- D$date
IsDate = function(x, format = NULL) {
  formatted = try(as.Date(x, format), silent = TRUE)
  is_date = as.character(formatted) == x & !is.na(formatted)  # valid and identical to input
  is_date[is.na(x)] = NA  # Insert NA for NA in x
  return(is_date)
}
dates_ix <- which(IsDate(dates, format = "%Y-%m-%d"))
D <- D[dates_ix,]
dates <- D$date

### filter data to start and end date
#dates <- D$date
u_dates <- unique(dates)
id_dates <- 1:length(u_dates)
between_dates <- subset(u_dates, (id_dates >= match(date_start, u_dates))&(id_dates <= match(date_end, u_dates)))
D <- D[D$date %in% between_dates, ]
sprintf("Timeframe of extracted mutation profiles %s to %s", unique(D$date)[1], unique(D$date)[length(unique(D$date))])

### Set empty lineage names to "nan":
#xi <- which(D$lineage == "")
#D$lineage[xi] <- "nan"

### filter mutations 
lineages <- sort(unique(D$lineage))
number_lineages <- length(lineages)
#print(D$date)

## check if lineages are in data time-horizon considered
#check_lin <- c('DV.6.1', 'DV.6.2', 'DV.7', 'DV.8', 'EG.5', 'EG.5.1.1', 'EG.5.2', 'EG.6.1', 'EG.7', 'FD.1.1', 'FE.1.1.1', 'FK.1', 'FK.1.2.1', 'FL.1.3', 'FL.1.5', 'FL.13', 'FL.2.3', 'FL.3.3', 'FL.6', 'FU.2', 'FY.1.2', 'FY.3', 'FY.4.1', 'GA.3', 'GE.1', 'GF.1', 'GG.1', 'GJ.1', 'GJ.1.1', 'GK.1', 'GN.1', 'GP.2', 'GR.1', 'XBB.1.16.11', 'XBB.1.16.6', 'XBB.1.33', 'XBB.1.34.1', 'XBB.1.42', 'XBB.1.47.1', 'XBB.1.5.28', 'XBB.1.5.68', 'XBB.1.5.70', 'XBB.1.5.72', 'XBB.1.5.73', 'XBB.1.5.75', 'XBB.1.5.86', 'XBB.2.3.11', 'XBB.2.3.8')
check_lin <- c("nan", "UNASSIGNED")

mutationprofiles_l <- list()
lineages_l <- c()
lineages_without_mutations <- c()
number_genomes_per_lineage <- c()
for (i in 1:number_lineages){
  D_lineagename <- lineages[i]
  Dlin <- D[which(D$lineage == D_lineagename),]
  D_N <- nrow(Dlin)
  aaprofile <- Dlin$aa_profile
  aaprofile[which(aaprofile=="")]<-NA
  #if ((length(which(!is.na(aaprofile))) > 0)&(D_lineagename != "")){ 
  if (length(which(!is.na(aaprofile))) > 0){ ## we need to consider all non-NaN aa_profiles
    aaprofile_l <- list()
    for (j in 1:D_N){
      aaprofile_l[j] <- strsplit(aaprofile[j]," ")
    }
    count <- table(unlist(lapply(aaprofile_l, unique)))
    count_df <- as.data.frame(count)
    count_df <- count_df[order(count_df$Freq, decreasing = TRUE),]
    N <- D_N*threshold
    D_mutationprofile <- count_df$Var1[which(count_df$Freq >= N)]
    mutationprofiles_l[i] <- list(D_mutationprofile)
    if (D_lineagename == ""){D_lineagename <- "nan"; lineages[i] <- "nan"} # renaming the relevant empty lineages to agree with downstream python scripts
    lineages_l = c(lineages_l, D_lineagename)
    number_genomes_per_lineage <- c(number_genomes_per_lineage, D_N)
    if (D_lineagename %in% check_lin){print(D_lineagename); print(D_mutationprofile)}
  }
  else{
    if (D_lineagename != ""){
      lineages_without_mutations <- c(lineages_without_mutations, D_lineagename) 
    }
  }
  #number_genomes_per_lineage <- c(number_genomes_per_lineage, D_N)
}
#number_lineages <- number_lineages-length(lineages_without_mutations)
number_genomes_per_lineage <- cbind(lineages_l, number_genomes_per_lineage)
write.csv(number_genomes_per_lineage, file=paste0(outputdir,"/",stringr::str_replace(outputfile_mutationprofile, ".csv","number_of_genomes_per_lineage.csv")), quote = FALSE, row.names = FALSE)
print(paste("Number of lineages in this dataset:",number_lineages))
rm(Dlin); rm(D_lineagename); rm(D_N);rm(aaprofile);rm(count);rm(count_df);rm(N);rm(D_mutationprofile)

## Grep for only spike ("S:") and RBD (RBDsites) mutations.
mutation_lists <- c()
for (i in 1:number_lineages){
  name <- lineages[i]
  m1 <- mutationprofiles_l[[i]]
  m1 <- m1[grep("S:",m1)]
  m1 <- gsub("S:","",m1)
  if (length(grep(":",m1) > 0)) {m1 <- m1[-grep(":",m1)] } #remove any deletions etc. keep only mutations
  m1 <- sort(m1)
  if (length(mutation_lists) == 0){
    mutation_lists <- cbind(name, paste(m1, collapse = "/"))
  }else{
    mutation_lists <- rbind(mutation_lists, 
                            cbind(name, paste(m1, collapse = "/")))
  }
}
colnames(mutation_lists) <- c("lineage", "mutated_sites_RBD")
mutation_lists <- as.data.frame(mutation_lists)
rm(m1)
write.csv(mutation_lists, file=paste0(outputdir,"/","mutation_lists.csv"), quote = FALSE, row.names = FALSE)

## remove empty entries in aa profile:
#ix_empty <- which(mutation_lists$mutated_sites_RBD == "")
#if (length(ix_empty) > 0){mutation_lists <- mutation_lists[-ix_empty,]}
#lineages <- mutation_lists$lineage
#number_lineages <- length(lineages)

## store as matrix:
spikemutations_list <- c()
for (i in 1:number_lineages){
  spikemutations_list <- c(spikemutations_list, strsplit(mutation_lists$mutated_sites_RBD[i],"/")[[1]])
}
spikemutations_list <- unique(spikemutations_list)

## Sort according to site:
spikemutations_list_site <- parse_number(spikemutations_list)
ix <- sort(spikemutations_list_site, index.return=TRUE)$ix
spikemutations_list <- spikemutations_list[ix]
rm(ix);rm(spikemutations_list_site)

MP <- as.data.frame(matrix(nrow = number_lineages, ncol = length(spikemutations_list)))
colnames(MP) <- spikemutations_list
rownames(MP) <- lineages
for (i in 1:number_lineages){
  m1 <- strsplit(mutation_lists$mutated_sites_RBD[i],"/")[[1]]
  x <- unlist(lapply(m1, function(x) which(colnames(MP) == x)))
  MP[lineages[i],x] <- 1
}
rm(m1);rm(x)

## Remove sites with empty entries only:
j <- c()
for (i in 1:ncol(MP)){
  if (length(which(!is.na(MP[,i]))) == 0){
    j <- c(j, i)
  }
}
if(length(j)>0){ 
  MP2 <- MP[,-j]
}else{
  MP2<- MP
}
MP2[is.na(MP2)]<-0
MP2 <- MP2[order(row.names(MP2)), ]
rm(MP)

write.csv(mutation_lists, file=paste0(outputdir,"/",outputfile_mutationlist), quote = FALSE, row.names = FALSE)
write.csv(MP2, file=paste0(outputdir,"/",outputfile_mutationprofile))





######################################### DATA PROCESSING 2: NTD & RBD
sites <- parse_number(colnames(MP2))
site_classification_rbd = ifelse(sites %in% RBD_positions, 1, 0)
site_classification_ntd = ifelse(sites %in% NTD_positions, 1, 0)
site_classification = site_classification_rbd + site_classification_ntd

mydf <- data.frame(row.names = colnames(MP2), category = site_classification)



######################################### DATA PROCESSING 3: Lineage Grouping

## Lineages are called identical if NTD and RBD regions are equal:
ix <- c(which(mydf==1))
MP3 <- MP2[,ix]
mydf3 <- as.data.frame(mydf[ix,]); colnames(mydf3) <- "category"; rownames(mydf3) <- rownames(mydf)[ix]
rm(ix)

## Some stats:
print(paste("Number of mutations in the spike protein for all given lineages:", ncol(MP2)))
print(paste("Number of mutations in the NTD region of the spike for all given lineages:", sum(site_classification_ntd)))
print(paste("Number of mutations in the RBD region of the spike for all given lineages:",sum(site_classification_rbd)))


## Spike-pseudogrouping: find duplicated lineages based on equal profile of NTD & RBD:
MP_l <- c()
for (i in 1:nrow(MP3)){
  MP_l <- c(MP_l, paste(colnames(MP3)[which(MP3[i,]==1)],collapse = "/"))
}
MP_l <- data.frame(cbind(rownames(MP3),MP_l))
colnames(MP_l) <- c("lineage","RBD_NTD_mutations")
write.csv(MP_l, file=paste0(outputdir,"/",outputfile_mutationprofile_mutations), row.names = FALSE, quote = FALSE)

psgroups <- c()
for (i in 1:nrow(MP_l)){
  #l <- MP_l$lineage[i]
  members <- MP_l$lineage[which(MP_l$RBD_NTD_mutations == MP_l$RBD_NTD_mutations[i])]
  psgroups <- c(psgroups, paste(members,collapse = "/"))
}
psgroups <- data.frame(cbind(MP_l$lineage,psgroups))
colnames(psgroups) <- c("lineage","group")

psgroups2 <- psgroups[!duplicated(psgroups$group),]
write.csv(psgroups2, file=paste0(outputdir,"/",outputfile_mutationprofile_groups), row.names = FALSE, quote = FALSE)


print(paste("Number of Spike-Pseudogroups found:", length(grep("/",psgroups2$group))))
print(paste("Number of individual lineages:",nrow(psgroups2) - length(grep("/",psgroups2$group))))


####################### Plotting:
# Limit the numer of lineages in the plotting: at least 100 genomes must be present for being plotted:
min_number_genomes_for_plotting <- 100
number_genomes_per_lineage <- data.frame(number_genomes_per_lineage)
colnames(number_genomes_per_lineage) <- c("lineage","N")
number_genomes_per_lineage$N <- as.numeric(number_genomes_per_lineage$N )
lineages_above_genome_threshold <- number_genomes_per_lineage$lineage[which(number_genomes_per_lineage$N > min_number_genomes_for_plotting)]

# Highlight Spike-Pseudogroups with "*" behind name to distinguish from individual lineages that do not belong to any group (or represent an own group, respectivly)
psgr <- c()
for (i in 1:nrow(psgroups2)){
  if (psgroups2$group[i] != psgroups2$lineage[i]){
    psgr <- c(psgr,psgroups2$lineage[i])
  }
}
MP3_red <- MP3[unlist(lapply(lineages_above_genome_threshold, function(x) which(rownames(MP3) == x))),]
MP4_unique <- MP3_red[!duplicated(MP3_red),]
x <- unlist(lapply(psgr, function(x) which(rownames(MP4_unique) == x)))
rownames(MP4_unique)[x] <- paste0(rownames(MP4_unique)[x],"*")
rm(x)

#remove empty mutational sites:
x <- which(apply(MP4_unique,2,sum)>0)
MP4_unique <- MP4_unique[,x]
mydf3 <- as.data.frame(mydf3[x,])
colnames(mydf3) <- "category"
rownames(mydf3) <- colnames(MP4_unique)
rm(x)
print(paste("Number of lineages / Spike-pseudogroups found with > ",min_number_genomes_for_plotting," genomes available in the dataset that are plotted as a heatmap:", nrow(MP4_unique)))

if (length(MP4_unique)>=2){
  pdf(paste0(outputdir,"/",outputfile_mutationprofile_plot), height = 15, width = 15)
  pheatmap(as.matrix(MP4_unique), main = paste("Spike Mutation Profile (NTD / RBD)"), col=c("white","red"),
           cluster_cols = F, cluster_rows = T,
           fontsize_col = 10,fontsize_row = 10,
           legend = FALSE,
           annotation_col = mydf3)
  dev.off()
}


## Zooming in for a predefined set of lineages:
lineages_tb <- data.frame(table(D$lineage))
# random selection:
# zoom_in_lineages <- sample(lineages_tb$Var1,3)
# take the top 3:
#zoom_in_lineages <- lineages_tb[order(lineages_tb$Freq, decreasing = TRUE)[1:3],1]
# take the top 3 of differentially mutated profiles of the RBD region:
MP3 <- MP2[,unlist(lapply(intersect(rownames(mydf3),colnames(MP2)),function(x) which(colnames(MP2) == x)))]
MP3_uniq <- unique(MP3)
lineags_tb2 <- subset(lineages_tb, Var1 %in% row.names(MP3_uniq))
zoom_in_lineages <- lineags_tb2[order(lineags_tb2$Freq, decreasing = TRUE)[1:3],1]

MP3_zoom <- MP3[unlist(lapply(zoom_in_lineages, function(x) which(rownames(MP3) == x))),]
#keep only those sites with at least one mutations among the lineages
MP3_zoom <- MP3_zoom[,which(apply(MP3_zoom,2,sum)>0)]
#keep only RBD and NTD regions:
#MP3_zoom <- MP3_zoom[,unlist(lapply(intersect(rownames(mydf3),colnames(MP3_zoom)),
#                                    function(x) which(colnames(MP2_zoom) == x)))]

if (length(MP3_zoom)>=2){
  pdf(paste0(outputdir,"/",stringr::str_replace(outputfile_mutationprofile_plot,".pdf","_zoom.pdf")),height = 3, width = 10)
  pheatmap(as.matrix(MP3_zoom), 
           main = paste("Zoom Mutation Profile of",paste(zoom_in_lineages,collapse = ",")), 
           col=c("white","red"),
           cluster_cols = F, cluster_rows = T,
           fontsize_col = 10,fontsize_row = 10,
           legend = FALSE,
           annotation_col = mydf3)
  dev.off()
}

