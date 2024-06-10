library(outbreakinfo)

lineages_of_interest <- "B.1.1.7"

############ Login to gisaid.org in browser to initialize the session in this script:
print("Please login to gisaid.org in browser to initialize the session in this script.")
outbreakinfo::authenticateUser()


mutations = getMutationsByLineage(pangolin_lineage=lineages_of_interest, frequency=0.75, logInfo = FALSE)

# Plot the mutations as a heatmap
#plotMutationHeatmap(mutations, title = "S-gene mutations in lineages")

# filter for Spike:
mutations_s <- subset(mutations, gene=="S")

# filter for RBD and NTD mutations:
#RBD_positions <- c(14:20, 140:158, 245:264, 331:531) 
mutations_s_rbd <- subset(mutations_s, codon_num%in% RBD_positions)

write.csv(mutations_s_rbd, "mutationprofile_outbreakinfo.csv", quote = FALSE, row.names = FALSE)

for (lin in lineages_of_interest){
   lin_mut <- subset(mutations_s_rbd, lineage==lin)
   lin_mut <- paste0(lin_mut$ref_aa,as.character(lin_mut$codon_num),lin_mut$alt_aa)
   write.table(lin_mut, paste0(dir_mutationprofiles_lineages,lin,"_mutations.txt"),quote = FALSE, row.names = FALSE,col.names = FALSE)
}
