rm(list = ls())
setwd("Desktop/percent_study/")
threshold = 3 #percent

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Call: Rscript Filter_lineages <lineage frequency data file> <filter threshold> <output directory>", call.=FALSE)
} else  {
  lineage_frequency_file <- args[1]
  threshold <- as.numeric(args[2])
  outputdir <- args[3]
}

## Outputfiles:
outputfile <- paste0(outputdir,"/filtered_lineages_")

# load the daily lineage frequencies in order to filter for variants who pass
# the threshold in the time frame
data_variant_percentage = read.csv(lineage_frequency_file)
data_variant_max = as.data.frame(apply(data_variant_percentage,2,max))

colnames(data_variant_max) = c("max", "variants")

dim(data_variant_max)

data_variant_max$variants = rownames(data_variant_max)
data_filtered = data_variant_max[data_variant_max$max > threshold, ]

variants = data_filtered$variants
length(variants)
colnames(variants) <- c("lineage")
### Save filtered variants
write.csv(variants, file=paste(outputfile, threshold,'_percent.csv'), row.names = FALSE, quote = FALSE)

