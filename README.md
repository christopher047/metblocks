
## metblocks
This package finds regions of variable methylation without group information. Package provides flexible parameters for finding small, large, densely packed, ... variably methylated regions called blocks. Package returns a list of blocks, imputed relative methylation matrix, and index to map blocks to matrix  

## installation
```R
#library(devtools)   
#library(GenomicRanges)   
#install_github("christopher047/metblocks")  
#library(metblocks) 
```

## Run metblocks on sample data chr18
library(GenomicRanges)  
library(metblocks)   

### return number of blocks, mean width in basepairs, number of cpgs, and cpg density per block
someResults <- function(blocks)  
        {  
        cat("number of blocks", length(blocks), "\n")  
        cat("mean block width in bp", mean(width(blocks)), "\n")  
        cat("mean number of cpgs per block", mean(blocks$L), "\n")  
        cat("cpgs per bp", round(mean(blocks$L)/mean(width(blocks)),3), "\n")  
        }   

### for parameter explanation see ?runChromosome

achr = "chr18"  
input_dir = file.path(system.file("extdata", package="metblocks"), "")   
min.seg = 20  
bwd = 300  
min.block = 5  
hclust = 0.3  
iqr_cutoff = 10  
nb = 5  
ncores = 16  

### paramaters used in paper 
cat("standard paramaters \n")   
chr18 <- runChromosome(achr, input_dir, min.seg, bwd, min.block, hclust, iqr_cutoff, nb, ncores)   
someResults(chr18$blocks)   

### Bandwidth (bwd) 

### smaller bandwidth gives smaller blocks more densely packed with cpgs 
cat("\n changing bandwidth to 150  \n")  
chr18 <- runChromosome(achr, input_dir, min.seg, bwd=150, min.block, hclust, iqr_cutoff, nb, ncores)   
someResults(chr18$blocks) 

### larger bandwidth gives larger blocks less densely packed with cpgs 
cat("\n changing bandwidth to 500\n")   
chr18 <- runChromosome(achr, input_dir, min.seg, bwd=500, min.block, hclust, iqr_cutoff, nb, ncores)   
someResults(chr18$blocks)

### Hclust (hclust)

### larger hclust gives more blocks that are less variable
cat("\n changing hclust to 0.4\n")   
chr18 <- runChromosome(achr, input_dir, min.seg, bwd, min.block, hclust=0.4, iqr_cutoff, nb,   ncores)   
someResults(chr18$blocks)  

### larger hclust gives fewer blocks that are more variable
cat("\n changing hclust to 0.2\n")  
chr18 <- runChromosome(achr, input_dir, min.seg, bwd, min.block, hclust=0.2, iqr_cutoff, nb, ncores)  
someResults(chr18$blocks)  

### Min.block (min.block)

### increasing min.block gives fewer blocks with more cpg sites
cat("\n changing min.block to 10\n")  
chr18 <- runChromosome(achr, input_dir, min.seg, bwd, min.block=10, hclust, iqr_cutoff, nb, ncores)  
someResults(chr18$blocks)  

### Min.seg (min.seg) 

### increasing min.seg will give larger segments with more cpg sites 
cat("\n changing min.seg to 40\n")  
chr18 <- runChromosome(achr, input_dir, min.seg=40, bwd, min.block, hclust, iqr_cutoff, nb, ncores)    
someResults(chr18$blocks)  


## Results 
blocks  <- chr18$blocks  #GRanges list of variable methylated regions or blocks         
index   <- chr18$index   #Index of CpG sites, blocks, and additional stats         
mat     <- chr18$m2      #The imputed matrix as data.frame       
segs    <- chr18$segs    #List of CpG sites in each segment    
raw_cov <- chr18$raw_cov #Number of reads per cite as a data.frame     

