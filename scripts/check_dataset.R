#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Call: Rscript check_dataset.R <covsonar data file> <country> <output file>", call.=FALSE)
} else  {
  input_datafile_covsonar <- args[1]
  input_country <- args[2]
  output_file <- args[3]
  #input_datafile_covsonar <- "gisaid_Australia_022022_112023.tsv"
  #input_country <- "Australia"
  #output_file <- "test.tsv"
}

D <- read.csv(input_datafile_covsonar, sep = "\t")
## Order by date:
D <- D[order(D$date),]

## Check country:
country <- unique(D$zip)
if (length(country) > 1){
  country_ix <- which(D$zip == input_country)
  if (length(country_ix)>0){
    D <- D[country_ix,]
  }
} 

## Check valid date format: "%Y-%m-%d"
IsDate = function(x, format = NULL) {
  formatted = try(as.Date(x, format), silent = TRUE)
  is_date = as.character(formatted) == x & !is.na(formatted)  # valid and identical to input
  is_date[is.na(x)] = NA  # Insert NA for NA in x
  return(is_date)
}
dates_ix <- which(IsDate(D$date, format = "%Y-%m-%d"))
D <- D[dates_ix,]

## Check for valid lineage entries:
lin_ix <- which(D$lineage == "UNASSIGNED")
lin_ix <- c(lin_ix, which(is.na(D$lineage)))
if(length(lin_ix)>0){
  D <- D[-lin_ix,]
}

## Check for valid mutation profile:
mutprofile_ix <- which(is.na(D$aa_profile))
mutprofile_ix <- c(mutprofile_ix, which(is.na(D$dna_profile)))
if(length(mutprofile_ix)>0){
  D <- D[-mutprofile_ix,]
}

write.table(D, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
