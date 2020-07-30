library(Biobase)
library(flowCore)
library(flowViz)
library(dplyr)
library(Seurat)

clr_normalized <- './analyzed/standard_workflow_besca2.0/citeseq/normalized_counts/clr/'
out_file <- './analyzed/standard_workflow_besca2.0/citeseq/citeseq.fcs'
mapping_file <-  './analyzed/standard_workflow_besca2.0/citeseq/mapping.csv'

#read input
h5ad <- Read10X(data.dir = clr_normalized, gene.column = 2, unique.features = TRUE)
dta <- t(as.matrix(h5ad))

#add column with barcode mapping
barcodes <- dimnames(dta)[[1]]
mapping <- data.frame(barcode_10X = barcodes, barcode_fcs = 1:length(barcodes))

#write mapping to file
write.table(x = mapping, file = mapping_file)

#add mapping to csv
dta <- cbind(dta, barcodes = mapping$barcode_fcs)

# you need to prepare some metadata
meta <- data.frame(name=dimnames(dta)[[2]],
                   desc=paste(dimnames(dta)[[2]])
)
meta$range <- apply(apply(dta,2,range),2,diff)
meta$minRange <- apply(dta,2,min)
meta$maxRange <- apply(dta,2,max)

head(meta)

# a flowFrame is the internal representation of a FCS file
ff <- new("flowFrame",
          exprs=dta,
          parameters=AnnotatedDataFrame(meta)
)

# now you can save it back to the filesystem
write.FCS(ff, out_file)


