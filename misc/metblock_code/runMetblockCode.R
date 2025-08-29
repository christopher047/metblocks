

### run metblocks ###
library(GenomicRanges)
library(metblocks)
set.seed(1024) 
input_dir <- file.path(system.file("extdata", package="metblocks"), "")
achr <- "chr18"
chr18 <- runChromosome(achr=achr, input_dir=input_dir)

### get results ###
blocks <- chr18$blocks #GRanges list of variable methylated regions or blocks
index <- chr18$index #Index of CpG sites, blocks, and additional stats
mat <- chr18$m2 #The imputed matrix as data.frame
segs <- chr18$segs #list of CpG sites in each segment
raw_cov <- chr18$raw_cov # number of reads per cite as a data.frame

### save results 
temp_block <- data.frame(blocks)
rownames(temp_block) <- names(blocks)
blocks <- temp_block
index <- data.frame(index)
rownames(index) <- paste0(index$seqnames, ".", index$start)

write.table(temp_block, "../results/raw_blocks_trim.tab", sep="\t")
write.table(index, "../results/raw_index_trim.tab", sep="\t")
write.table(mat, "../results/raw_mat_trim.tab", sep="\t")


