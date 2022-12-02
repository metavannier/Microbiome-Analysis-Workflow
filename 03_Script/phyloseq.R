# #################################################################
# 
#        Construct the phyloseq object
#          
# #################################################################

library(phyloseq)
library(qiime2R)
library(dplyr)
library(readr)

# Input files
filtertable <- snakemake@input[["filtertable"]]
rooted_tree <- snakemake@input[["rooted_tree"]]
sklearn <- snakemake@input[["sklearn"]]
metadata <- snakemake@input[["metadata"]]
output_phylseq <- snakemake@output[["output_phylseq"]]

physeq <- qza_to_phyloseq(
  features=filtertable,
  tree=rooted_tree,
  sklearn,
  metadata = metadata
)
physeq
#write.table(physeq, file = output_phylseq, sep="\t", quote = FALSE,row.names = FALSE, col.names = FALSE)


