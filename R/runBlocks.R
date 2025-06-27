
#### create object and run ###
source("./R/block_methods.R")
source("./R/blockClassDefinition.R")
require(GenomicRanges) 

#' runChromosome
#'
#' Main function for finding blocks on a chromosome
#' @param achr a name of a chromosome for example chr18   
#' @param input_dir directory containing chr18_tot_cov.tab and chr18_tot_met.tab files  
#' @param min.seg is the minimum number of CpG sites needed to keep the segment 
#' @param bwd or distance decay bandwith parameter 
#' @param min.block minimum number of CpG site per block 
#' @param hclust the cutoff used in hclust to divide into blocks   
#' @param iqr_cutoff the range/iqr to remove blocks with 1-2 outliers  
#' @param nb number of neighbours used by KNN to impute missing values 
#' @param ncores number of cores to use 
#' @return a list of results including blocks with statistics  
#' @export

runChromosome <- function(achr="chr18", input_dir="./inst/extdata/", min.seg=20, bwd=300, min.block=5, hclust=0.3, iqr_cutoff=10, nb=5, ncores=16) 
	{
	test <- createMethObject() 
	test <- test@readFiles(test, input_dir=input_dir, achr=achr) 
	test <- test@createIndex(test) 
	test <- test@getSegs(test, min.seg=min.seg) 
	test <- test@imputeKNNbySegs(test, ncores=ncores, nb=nb) 
	test <- test@getBlocks(test, bwd=bwd, min.block=min.block, hclust=hclust, ncores=ncores)
	test <- test@filterIQR(test, iqr_cutoff = iqr_cutoff) 
	test <- test@calcBlockCoverage(test, min.block=min.block, ncores=ncores) 
	res  <- test@getReturnData(test)
	if(!all(!is.na(test@index))){return(c())}
        return(res) 
	}



#' test_metblocks
#'
#' Main function runs shows an example on chr18
#' @return None
#' @export


test_metblocks <- function() 
	{
	#set.seed(1024)
	#input_dir <- paste0(system.file("extdata", package="metblocks"), "/")  
	#chr18         <- runChromosome(achr="chr18", input_dir=input_dir)
	#blocks        <- chr18$blocks
	#index         <- chr18$index
	#imp_mat       <- chr18$m2
	#segs          <- getSegRanges(chr18$segs)
	#cov           <- chr18$raw_cov

	#### write blocks ####
	#temp_block           <- data.frame(blocks)
	#rownames(temp_block) <- names(blocks)
	#write.table(temp_block, "../results/raw_blocks_trim.tab", sep="\t")

	#### write index ####
	#index           <- data.frame(index)
	#rownames(index) <- paste0(index$seqnames, ".", index$start)
	#write.table(index, "../results/raw_index_trim.tab", sep="\t")

	### write imputed matrix ###
	#mat <- chr18$m2
	#write.table(mat, "../results/raw_mat_trim.tab", sep="\t")

	### write segments #####
	#write.table(data.frame(segs), "../results/segs_trim.tab", sep="\t")
	}



