# metblocks

# installation
library(devtools)   
library(GenomicRanges)   
install_github("christopher047/metblocks")  
library(metblocks)    

# Run metblocks on sample data chr18
library(GenomicRanges)  
library(metblocks)   
input_dir <- file.path(system.file("extdata", package="metblocks"), "")   
achr <- "chr18"  
chr18 <- runChromosome(achr=achr, input_dir=input_dir)   

# Results 
blocks <- chr18$blocks #GRanges list of variable methylated regions or blocks       
index  <- chr18$index #Index of CpG sites, blocks, and additional stats       
mat    <- chr18$m2 #The imputed matrix as data.frame     
segs   <- chr18$segs #list of CpG sites in each segment     

# Saving Results

temp_block           <- data.frame(blocks)  
rownames(temp_block) <- names(blocks)  
blocks               <- temp_block  
index                <- data.frame(index)  
rownames(index)      <- paste0(index$seqnames, ".", index$start)  

