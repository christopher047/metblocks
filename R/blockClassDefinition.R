
#' setClass
#'
#' This function defines the meth class attributes 
#' @return an meth class definition
#' @export  

### this is just for swapping out methods under development ####
setClass("meth", slots=list(achr    = "character",  #the chromosome being processed
                            imp_mat ="data.frame",  #the imputed matrix 
                            raw_cov = "data.frame", #the raw coverage 
                            raw_meth = "data.frame",#the raw methylation
                            segs    ="list",        #the segments 
                            index   ="GRanges",     #GRange index of each site
                            blocks  ="GRanges",     #GRange index of each block
                            readFiles ="function",  #reads the files 
                            createIndex = "function", #create index  
                            getSegs   = "function", #divided chromosomes into segments 
                            imputeKNNbySegs = "function",#imputes the imp_matrix by segment, imp_mat is raw_meth/raw_counts, raw_count =0 gives NA 
                            getBlocks       = "function", #get blocks by variation and distance
                            filterIQR       = "function", #filters blocks with one or two outliers 
                            calcBlockImputation = "function", #calculates the percent impution per block
                            calcBlockCoverage   = "function", #calculates L(number of cpgs), block_sum(number total reads), block_mean(mean relative methylation)
                            getReturnData       = "function"))# returns blocks and imp_mat of blocks 

#' createMethObject
#'
#' This constructs a meth class object with attributes
#' @return an meth class object
#' @export  

createMethObject <- function()
        {
        temp <- new("meth")
        temp@readFiles           <- ReadRawFiles
        temp@createIndex         <- createIndex
        temp@getSegs             <- getSegs
        temp@imputeKNNbySegs     <- imputeKNNbySegs
        temp@getBlocks           <- getBlocks
        temp@filterIQR           <- filterIQR
        temp@calcBlockCoverage   <- calcBlockCoverage
        temp@getReturnData       <- getReturnData
        return(temp)
        }

