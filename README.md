# metblocks

# installation
library(devtools)   
library(GenomicRanges)   
install_github("christopher047/metblocks")  
library(metblocks)    

# Run metblocks on sample data chr18
input_dir <- file.path(system.file("extdata", package="metblocks"), "")   
achr <- "chr18"  
chr18 <- runChromosome(achr=achr, input_dir=input_dir)   

# Results 
