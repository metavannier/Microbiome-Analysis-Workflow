# #################################################################
# 
#       Table of the relative abundance after taxa colapse, 
#       feature filtering per group.
#       Format is for the input of LEfSe tools
#              
# #################################################################

library(readr)
library(formattable)
# library(tibble)

# Output table with LEfSe format
table_abond_lefse=snakemake@output[["table_abond_lefse"]]

# Relative abudance table to transform as LEfSe format
table_abond_tsv <- read.table(snakemake@input[["table_abond_tsv"]], header=FALSE, sep = "\t")

table_abond_tsv[,1] <- gsub('\"','',table_abond_tsv[,1])

# Add marker column
table_abond_colname <- read.csv(snakemake@input[["table_abond_tsv"]], header=FALSE, sep = "\t")
table_abond_colname <- table_abond_colname[-1,-ncol(table_abond_colname)]
table_abond_colname[1,1] <- "index"

# Match with metadata to have the group corresponding to the sample name
table_abond_colname2 <- table_abond_colname[,-1]
metadata <- read.delim(snakemake@params[["metadata"]], header=TRUE, sep = "\t")
# Sample in the same order than the original abundance table
group_idx <-  match(table_abond_colname2[1,],metadata$sample.id)
group_experiment  <- metadata[group_idx,ncol(metadata)]

# Add the rowname of the group experiment row
group_experiment <- append("day", group_experiment)

# ADD the index and day row to the original abundance datatable
index_day <- rbind(table_abond_colname[1,],group_experiment)
table_abond_lefse_final <- rbind(index_day,table_abond_colname)
table_abond_lefse_final <- table_abond_lefse_final[-3,]
table_abond_lefse_final[,1] <- gsub(';','|',table_abond_lefse_final[,1])

write.table(table_abond_lefse_final, file = table_abond_lefse, sep="\t", quote = FALSE,row.names = FALSE, col.names = FALSE)
