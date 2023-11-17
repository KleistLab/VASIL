rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Call: Rscript Filter_lineages <lineage frequency data file> <filter threshold> <output directory>", call.=FALSE)
} else  {
  lineage_frequency_file <- args[1]
  threshold <- as.numeric(args[2])
  outputfile <- args[3]
}

# load the daily lineage frequencies in order to filter for variants who pass
# the threshold in the time frame
data_variant_percentage = read.csv(lineage_frequency_file)
drop <- c("X", "date", "week_num")
data_variant_percentage <- data_variant_percentage[, !(names(data_variant_percentage)%in%drop)]

#data_variant_max <- as.data.frame(apply(data_variant_percentage,2,max)) # does not work on Monterey
#colnames(data_variant_max) <- c("max", "variants") # does not work on Monterey
#data_variant_max$variants = rownames(data_variant_max)
#data_filtered = data_variant_max[data_variant_max$max > threshold, ]
#dim(data_variant_max)
#print(variants)
#variants = data_filtered$variants
#length(variants)

all_lineages <- colnames(data_variant_percentage)
variants <- c()
max_list <- c()
for (i in 1:length(all_lineages)){
  max_per <- max(data_variant_percentage[, all_lineages[i]])
  if (max_per > threshold)
    {
     variants <- append(variants, all_lineages[i])
     max_list <- append(max_list, max_per)
    }
}

variants_df <- data.frame("lineage"=variants, "max"=max_list)
### Save filtered variants
write.csv(variants_df, file = outputfile, row.names = FALSE, quote = FALSE)

