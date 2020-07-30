#read FCS file
library(Biobase)
library(flowCore)
library(flowViz)
library(dplyr)

path <- './celltypes/fcs'
outdir <- './celltypes/barcodes/'
prefix <- 'export_null_'
mapping <- "analyzed/standard_workflow_besca2.0/citeseq/mapping.csv"

#read mapping file
mapping <- read.csv(mapping, sep = " ")
barcodes <- mapping$barcode_10X
names(barcodes) <- mapping$barcode_fcs

#define lookup function
lookup <- function(x){
  value <- unname(barcodes[x])
  return(value)
}

for (file in list.files(path)) {
  print(file)
  name <- basename(file)
  name <- gsub(".fcs", ".csv", name)
  name <- gsub(prefix, "", name)
  
  #import FCS file
  samp <- read.FCS(paste0(path, '/', file), ignore.text.offset = TRUE)
  
  #convert to csv file
  samp.df <- as.data.frame(exprs(samp))
  
  #map ids to barcodes
  samp.df$barcodes_10X <- lapply(list(samp.df$barcodes), lookup)

  #extract barcodes and export to file
  write.table(as.data.frame(samp.df$barcodes_10X), 
              paste0(outdir, '/', name),
              row.names=FALSE,
              col.names = FALSE,
              sep = ',',
              append = FALSE)
}
