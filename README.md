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
blocks <- getResultBlocks(chr18) #GRanges list of variable methylated regions or blocks     
index  <- getResultIndices(chr18) #Index of CpG sites, blocks, and additional stats     
mat    <- getResultMatrices(chr18) #The imputed matrix    
segs   <- getResultSegs(chr18) #Grange list of segments   

# Saving Results

temp_block           <- data.frame(blocks)  
rownames(temp_block) <- names(blocks)  
blocks               <- temp_block  
index                <- data.frame(index)  
rownames(index)      <- paste0(index$seqnames, ".", index$start)  
segs                 <- data.frame(apply(data.frame(segs), 2, as.character))  

