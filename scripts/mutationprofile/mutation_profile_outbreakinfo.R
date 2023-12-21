library(outbreakinfo)
library(stringr)
library(tidyr)
library(dplyr)
library(readr)


outputdir <- "outbreakinfo/"
output_prefix_for_file <- "outbreakinfo"
dms_per_ab_per_site_file <- "dms_per_ab_per_site.csv"
dir.create(outputdir, showWarnings = FALSE)



############ Login to gisaid.org in browser to initialize the session in this script:
print("Please login to gisaid.org in browser to initialize the session in this script.")
outbreakinfo::authenticateUser()



# Define outputfiles:
outputfile_mutations_outbreakinfo <- paste0(outputdir,"/",output_prefix_for_file,"_downloaded_mutations_from_outbreakinfo.csv")
outputfile_mutationlist <- paste0(outputdir,"/",output_prefix_for_file,"_mutations_spike_lists.csv")
outputfile_mutationprofile <- paste0(outputdir,"/",output_prefix_for_file,"_mutations_spike.csv")
outputfile_mutationprofile_mutations <- paste0(outputdir,"/",output_prefix_for_file,"_RBD_NTD_mutations.csv")
outputfile_mutationprofile_groups <- paste0(outputdir,"/",output_prefix_for_file,"_RBD_NTD_pseudogroups.csv")
outputfile_stats <- paste0(outputdir,"/",output_prefix_for_file,"_stats.txt")

######## CONSTANT DEFINTIONS:
## spike protein positions (NTD and RBD):
dms_per_ab_per_site <- read.csv(dms_per_ab_per_site_file, sep = ",")
RBD_positions <- sort(unique(dms_per_ab_per_site$site))
RBD_position_start = min(RBD_positions)
RBD_position_end = max(RBD_positions)
NTD_positions <- c(14:20, 140:158, 245:264)
NTD_position_start = min(NTD_positions)
NTD_position_end = max(NTD_positions)




############  Download all current pango lineages from github
## https://github.com/cov-lineages/pango-designation/blob/master/lineage_notes.txt
url <- "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt"
lineages_pango_file <- paste(outputdir, "lineage_notes.txt", sep = "/") #"lineage_notes.txt"
## If old pango file exists, remove it and replace it with new one from the download:
if (file.exists(lineages_pango_file)) {
  file.remove(lineages_pango_file)
}
print("Downloading current pango lineages from github.")
download.file(url, lineages_pango_file, mode = "wb")



############ Read in all pango lineages and get mutation info from outbreak.info:
lineages_pango <- read.table(lineages_pango_file,sep = "\t",header = TRUE)
lineages_pango <- unique(lineages_pango[,1])

lineages_of_interest <- lineages_pango  #c("BA.2", "BA.2.12.1", "BA.4", "BA.5")
number_lineages <- length(lineages_of_interest)

# Get all mutations in the lineages of interest with at least 75% prevalent in one of the lineages.
print("Retrieving mutations by lineage from outbreak.info. This might take some minutes...")
mutations = getMutationsByLineage(pangolin_lineage=lineages_of_interest, frequency=0.75, logInfo = FALSE)
write.csv(mutations, file = outputfile_mutations_outbreakinfo, quote = FALSE, row.names = FALSE)


############ Plotting 
# Plot the mutations as a heatmap from outbreak.info
# plotMutationHeatmap(mutations, title = "S-gene mutations in lineages")



############ Extract Spike protein and generate mutation profiles for each lineage
## Parts of the code (editing the outbreak.info file) are extracted from https://github.com/outbreak-info/R-outbreak-info/blob/main/R/plotMutationHeatmap.R
if(!is.null(mutations) && nrow(mutations) != 0){
  mutations = mutations %>% filter(gene == "S")
  mutations = mutations %>%
    rowwise() %>%
    mutate(s_mutation = toupper(str_split(mutation, ":")[[1]][2])) %>%
    arrange(codon_num)
  
  # create empty grid
  s_mutation = mutations %>% pull(s_mutation) %>% unique()
  lineage = mutations %>% pull(lineage) %>% unique()
  blank = crossing(lineage, s_mutation)
  
  # refactor the mutations to sort them
  blank$s_mutation = factor(blank$s_mutation, levels = s_mutation)
  mutations$s_mutation = factor(mutations$s_mutation, levels = s_mutation)
  
  mutations2 <- mutations[,c("s_mutation","lineage","prevalence")]
  #Remove deletions:
  mutations2 <- mutations2[-grep("DEL",mutations2$s_mutation),]
  #Extract sites:
  mutations2$site <- parse_number(as.character(mutations2$s_mutation))
  
  ## Store as list:
  mutation_lists <- c()
  for (i in 1:number_lineages){
    xi <- which(mutations$lineage == lineages_of_interest[i])
    m1 <- mutations2$s_mutation[xi]
    m1 <- m1[!is.na(m1)]
    m1 <- unique(sort(m1))
    if (length(mutation_lists) == 0){
      mutation_lists <- cbind(lineages_of_interest[i], paste(m1, collapse = "/"))
    }else{
      mutation_lists <- rbind(mutation_lists, 
                              cbind(lineages_of_interest[i], paste(m1, collapse = "/")))
    }
  }
  colnames(mutation_lists) <- c("lineage", "mutated_sites_spike")
  mutation_lists <- as.data.frame(mutation_lists)
  
  ## store as matrix:
  spikemutations_list <- c()
  for (i in 1:number_lineages){
    spikemutations_list <- c(spikemutations_list, strsplit(mutation_lists$mutated_sites_spike[i],"/")[[1]])
  }
  spikemutations_list <- unique(spikemutations_list)
  
  ## Sort according to site:
  spikemutations_list_site <- parse_number(spikemutations_list)
  ix <- sort(spikemutations_list_site, index.return=TRUE)$ix
  spikemutations_list <- spikemutations_list[ix]
  rm(ix);rm(spikemutations_list_site)
  
  MP <- as.data.frame(matrix(nrow = number_lineages, ncol = length(spikemutations_list)))
  colnames(MP) <- spikemutations_list
  rownames(MP) <- lineages_of_interest
  for (i in 1:number_lineages){
    m1 <- strsplit(mutation_lists$mutated_sites_spike[i],"/")[[1]]
    x <- unlist(lapply(m1, function(x) which(colnames(MP) == x)))
    MP[lineages_of_interest[i],x] <- 1
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
  
  write.csv(mutation_lists, file=outputfile_mutationlist, quote = FALSE, row.names = FALSE)
  write.csv(MP2, file=outputfile_mutationprofile)
  
  
  ##############  Extracting NTD & RBD
  sites <- parse_number(colnames(MP2))
  site_classification_rbd = ifelse(sites %in% RBD_positions, 1, 0)
  site_classification_ntd = ifelse(sites %in% NTD_positions, 1, 0)
  site_classification = site_classification_rbd + site_classification_ntd
  
  mydf <- data.frame(row.names = colnames(MP2), category = site_classification)
  
  ############## Pseudo grouping
  ## Lineages are called identical if NTD and RBD regions are equal:
  ix <- c(which(mydf==1))
  MP3 <- MP2[,ix]
  mydf3 <- as.data.frame(mydf[ix,]); colnames(mydf3) <- "category"; rownames(mydf3) <- rownames(mydf)[ix]
  rm(ix)
  
  ## Some stats:
  st1 <- paste("Number of mutations in the spike protein for all given lineages:", ncol(MP2))
  print(st1); cat(st1, file = outputfile_stats, sep = "\n")
  st2 <- paste("Number of mutations in the NTD region of the spike for all given lineages:", sum(site_classification_ntd))
  print(st2); cat(st2, file = outputfile_stats, sep = "\n", append = TRUE)
  st3 <- paste("Number of mutations in the RBD region of the spike for all given lineages:",sum(site_classification_rbd))
  print(st3); cat(st3, file = outputfile_stats, sep = "\n", append = TRUE)
  
  
  ## Spike-pseudogrouping: find duplicated lineages based on equal profile of NTD & RBD:
  MP_l <- c()
  for (i in 1:nrow(MP3)){
    MP_l <- c(MP_l, paste(colnames(MP3)[which(MP3[i,]==1)],collapse = "/"))
  }
  MP_l <- data.frame(cbind(rownames(MP3),MP_l))
  colnames(MP_l) <- c("lineage","RBD_NTD_mutations")
  write.csv(MP_l, file=outputfile_mutationprofile_mutations, row.names = FALSE, quote = FALSE)
  
  psgroups <- c()
  for (i in 1:nrow(MP_l)){
    members <- MP_l$lineage[which(MP_l$RBD_NTD_mutations == MP_l$RBD_NTD_mutations[i])]
    psgroups <- c(psgroups, paste(members,collapse = "/"))
  }
  psgroups <- data.frame(cbind(MP_l$lineage,psgroups))
  colnames(psgroups) <- c("lineage","group")
  
  psgroups2 <- psgroups[!duplicated(psgroups$group),]
  write.csv(psgroups2, file=outputfile_mutationprofile_groups, row.names = FALSE, quote = FALSE)
  
  st4 <- paste("Number of Spike-Pseudogroups found:", length(grep("/",psgroups2$group)))
  print(st4); cat(st4, file = outputfile_stats, sep = "\n", append = TRUE)
  st5 <- paste("Number of individual lineages:",nrow(psgroups2) - length(grep("/",psgroups2$group)))
  print(st5); cat(st5, file = outputfile_stats, sep = "\n", append = TRUE)
} 
