rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Call: Rscript Filter_lineages <covsonar data file> <lineage frequency data file> <filter threshold> <output file filtered variants> <output file for filtered covsonar data>", call.=FALSE)
} else  {
  covsonar_data <- args[1]
  lineage_frequency_file <- args[2]
  threshold <- as.numeric(args[3])
  Var_outputfile <- args[4]
  Covsonar_outputfile <- args[5]
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
write.csv(variants_df, file = Var_outputfile, row.names = FALSE, quote = FALSE)

# load the stichproben file and delete all variants, that are below the threshold

data_stichprobe = read.csv(covsonar_data, sep = "\t")
head(data_stichprobe)

data_stichprobe_filtered = data_stichprobe[data_stichprobe$lineage %in% variants,]

date_max = max(data_stichprobe_filtered$date)
date_min = min(data_stichprobe_filtered$date)

date_min
date_max

# order of dates might be corrupted. If dates are not in order, VASIL crashes
data_stichprobe_filtered <- data_stichprobe_filtered[order(data_stichprobe_filtered$date),]

write.table(data_stichprobe_filtered, file=Covsonar_outputfile, quote=FALSE, sep='\t')