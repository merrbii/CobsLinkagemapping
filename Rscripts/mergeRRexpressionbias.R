# This script merges dataset published along the Schrader et al. 2017 manuscript 
# With file containing gene IDs in Cobs3.1 + gwRR

# load library
library(tidyverse)

# load Schrader data
schrader <- read_delim("RNAseq/msw240_Supp/SupplementaryTable.csv", col_names = T, delim = ",")

# load file with Cobs3.1 gene ids and corresponding Cobs1.4 ids
mygenes <- read_delim("RNAseq/Cobs3.1.geneID.Cobs1.4.geneID.txt", col_names = F, delim = "\t")

colnames(mygenes) <- c("cobs31","chr", "start", "end","geneID")

merged <- merge(mygenes, schrader, by = "geneID", all.x = T)

## Add recombination rates:
# load rr generated previously see Rscripts/gwRecombinationRate.R
rr <- read_delim("windowedRR/finalRecombinationRateEstimates.250kbwindows.txt", col_names = T)
rr <- rr[,c(1,11)]
colnames(rr) <- c("window_name", "rr")


# load gene window names: a bed file containing 5 columns: 	 <geneID(Cobs3.1)> <chr> <start> <end> <window_name>
genewid <- read_delim("RNAseq/genesCobs3.1.matching.ExpressedgenesCobs1.4.RR.tmp.txt", col_names = F)
colnames(genewid) <- c("cobs31", "chr", "start", "end", "window_name")

# combine datasets
merged2 <- merge(genewid, rr, by = "window_name")
merged3 <- merge(merged, merged2, by = c("cobs31", "chr", "start", "end"))

# remove entries with no expression data:
merged3 <- merged3[!is.na(silver$`Expr(Ave)`),]

# remove duplicate genes: i.e. genes split in the new assembly
n_occur1 <- data.frame(table(merged3$geneID))
merged3 <- merged3[merged3$geneID %in% n_occur1$Var1[n_occur1$Freq < 2],]

data <- merged3

# store table for correlation analysis:
write.table(data,"RNAseq/finals/genesCobs3.1.matching.ExpressedgenesCobs1.4.Schraderetal.withRR.txt", quote = F, sep = "\t", row.names = F, col.names = T)
