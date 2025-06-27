#' Load coverage and methylation files 
#'
#' This function loads two files (total coverage and total methylated). 
#' The coverage tot_cov.tab is the number of reads at each CpG site 
#' The methylation file is the number of methylated reads at each CpG site 
#' The naming convention for the coverage file is  chr??_tot_cov.tab
#' The naming convention for the methylation file is chr??_tot_meth.tab 
#' Where ?? is the chrosome (chr17) and input directory(./data) 
#  This point to files ./data/chr17_tot_cov.tab and ./data/chr17_tot_meth.tab
#' Files must be tab delimited 
#' The rownames must be in the following format chr.pos for example chr17.10002
#' The colnames contains are the sample identifiers.
#' Duplicate rownames or colnames are not allowed
#' Chromsome is stored in object@chr 
#' Total coverage is stored in object@raw_cov
#' Total methylated is stored in object@raw_meth
#' @param object a meth class object 
#' @param input_dir the directory of the count and methylation files 
#' @param achr is the name of the chromosome for example chr17
#' @return A meth class object 
#' @export 	

ReadRawFiles <- function(object, input_dir="./", achr="chr17")
                {
		#### reads raw coverage and raw methylation files ####	
		require(data.table) 	
		my_read <- function(filename) 
			{
			ff <- suppressWarnings(data.frame(fread(filename,sep="\t"))) 
                	rownames(ff) <- ff$V1
                	ff <- ff[,-1]
                	return(ff)
			}
		cat(paste("running", achr,"\n"))  
		object@achr	    <- achr
		object@raw_cov      <- my_read(paste0(input_dir, achr, "_tot_cov.tab"))  
        	object@raw_meth     <- my_read(paste0(input_dir, achr, "_tot_meth.tab")) 
		return(object) 
                }

#' createIndex  
#'
#' This function creates converts rownames chr17.1002 to a Granges object  
#' Ensures that the rownames are equal and sorted 
#' Generates rowsums for each CpG site  
#' Generated rowmeans for each CpG site
#' Results are stored in object 
#' The naming convention for the coverage file is  chr??_tot_cov.tab
#' The naming convention for the methylation file is chr??_tot_meth.tab 
#' Files are table delimited 
#' The rownames must be in the following format chr.pos for example chr17.10002
#' The colnames contains are the sample identifiers.
#' Duplicate rownames are not allowed
#' Modified total coverage is stored in object@raw_reads
#' Modified Total methylated is stored in object@raw_meth
#' Granges index is stored in object@index 
#' @param object a meth class object 
#' @param input_dir the directory of the count and methylation count files 
#' @param achr is the name of the chromosome for example chr17
#' @return A meth class object
#' @export   



createIndex  <- function(object) 
		{
		######## paramaters #####	
		tot_cov  <- object@raw_cov
		tot_meth <- object@raw_meth 
		#### keep coverage in index ####
		cat("creating index and row coverage\n") 
		chr            <- gsub("\\..*$","", rownames(tot_cov))
                pos            <- gsub("^.*\\.","", rownames(tot_cov))
                index          <- GRanges(data.frame("chr"=chr, "start"=pos, "end"=pos))
		names(index)   <- paste0(chr, ".", pos) 
        	########### create index ########
        	index            <- index[rownames(tot_cov)]
		index            <- sort(index)
		index$site_reads <- rowSums(tot_cov[names(index),], na.rm=T)  
        	index$site_mean  <- rowMeans(tot_cov[names(index),], na.rm=T) 	
		############ make sure they are in the same order ##########
        	object@raw_cov  <- tot_cov[names(index),]
        	object@raw_meth <- tot_meth[names(index),]
		object@index    <- index 
		return(object) 
		}

#' getSegs
#'
#' This function divides the chromome into segments based on distance 
#' This function filters out segment less than min.seg parameter 
#' Distance paramater is the closest distance to 3000 bp 
#' Segments are stored in the as GRanges in object@segs 
#' @param object of class meth  
#' @param min.seg is the minimum number of CpG sites needed to keep the segment 
#' @return an object of meth class 
#' @export  


getSegs <- function(object, min.seg=5)
        {
	#### paramaters ####
	index     <- object@index
	#### get closest cutoff to 3000 ####
        test  <- diff(start(index))
        keep <- as.numeric(gsub("%", "", (names(which.max(which(quantile(test, seq(0,1,0.01)) < 3000))))))/100
        cutoff <- as.integer(quantile(test, keep))
        cat(paste("segment cutoff is", cutoff, "\n"))
        cat(paste("quantile is", keep, "\n"))
        #setup segmentation list
        test <- c(0, test)
        test[test >  cutoff] <- 0
        test <- setNames(test, as.character(names(index)))
        segs <- vector(mode="list", sum(test==0))
        #segment change to GRanges
        j = 0
        #index <- c()
        for(i in 1:length(test))
                {
                if (test[i] == 0){j <- j+1}
                segs[[j]] <- c(segs[[j]], names(test)[i])
                }
        #keep only groups of min.block or more
        segs <- lapply(segs, na.omit)
        segs <- segs[which(as.numeric(lapply(segs, length)) >=min.seg)]
	object@segs <- segs 
	return(object) 
        }

#' getGrangeIndex 
#'
#' Helper function for converting row.names ie chr17.10001 to Granges  
#' @param xx a vector of characters ie c("chr17.10001, chr17.10002, ...) 
#' @param achr a chromosome name "chr17"
#' @return a sorted Rrange
#' @export   

getGRangeIndex <- function(xx, achr)
        {
	### this takes chromosome, smallest, largest to create a GRange	
        st    <- lapply(xx, function(i){min(as.integer(gsub(".*\\.", "", i)))})
        ed    <- lapply(xx, function(i){max(as.integer(gsub(".*\\.", "", i)))})
        gg    <- GRanges(data.frame("chr"=achr, start=as.integer(st), end=as.integer(ed)))
        return(sort(gg))
	### uncomment to remove overlapping blocks caused by hclust ###
	#return(reduce(sort(gg)))  
        }

#' getSegRanges 
#'
#' Helper function for converting segments to Granges
#' Adds number of methylation events and length 
#' expects a list of segments 
#' @param segs  
#' @return a table of segment with methylation events and length
#' @export 


getSegRanges <- function(segs) 
	{	
	### this function converts segs to GRanges adds number of M sites as L	
	achr   <- gsub("\\..*", "", segs[[1]][1])
	temp   <- getGRangeIndex(segs, achr) 
        temp$L <- lapply(segs, length)
	temp   <- apply(data.frame(temp), 2, as.character)
	return(temp)
	}
#' getSegIndex
#'
#' Helper function creates list of segment names and corresponding GRange
#' @param segs a list of segments   
#' @param achr a character chromsome name ie chr17  
#' @return list with segment names as keys and GRange as values
#' @export

createSegIndex <- function(segs, achr) 
	{
	### this names range by chromosome and start position	
	names(segs)   <- paste0(achr, ".", 1:length(segs))
        sindex        <- getGRangeIndex(segs, achr)
        names(sindex) <- paste0(achr, ".", 1:length(segs))
	return(sindex) 
	}	

#' getimputeKNNbySegs
#'
#' This function imputes missing values by KNN per segment
#' Stores the imputed relative methylation values (raw_meth/raw_count) in object@imp_mat
#' @param object an object of meth class
#' @param ncores the number of cores to use 
#' @param nb the number of neighhours used in KNN
#' @return an object of meth class
#' @export

imputeKNNbySegs <- function(object, ncores=16, nb=5) 
	{
	cat("imputing by segment \n") 	
	require(impute)
	require(parallel) 
	segs <- object@segs 	
	temp <- object@raw_meth/object@raw_cov
        x1       <- mclapply(1:length(segs), function(i){impute.knn(t(temp[segs[[i]],]), k=nb)$data}, mc.cores=ncores)
        b0       <- data.frame(t(do.call("cbind", x1)))
	object@imp_mat <- b0
	return(object) 
	}

#' methylDist
#'
#' Calculates a correlation matrix of segment relative methylation d1=cor(seg)
#' Calculates a distance matrix by segment position and weighted by bandwidth
#' d2 = exp(-as.matrix(dist(pos)^2)/(2*bwd^2))
#' Returns d, where d=1 - d1*d2  
#' @param X a segment 
#' @param pos numeric postion of methylation sites on segment 
#' @param bwd bandwidth
#' @return a matrix generated from covariation, and weighted by distance
#' @export


methylDist <- function(X,pos,bwd = 300)
        {
        d1 = cor(t(X))
        d1[is.na(d1)] <- 1
        d2 <- exp(-as.matrix(dist(pos)^2)/(2*bwd^2))
        d <- 1-(d1*d2)
        return(d)
        }

#' endreCluster
#'
#' Divides segments into blocks using methyDist and hclust
#' Returns list of blocks per segment 
#' Filters out blocks with less than min.block CpG sites 
#' @param bb a segment 
#' @param bwd bandwidth
#' @param min.block minumum number CpGs sites per block
#' @param hclust hclust cutoff for dividing blocks
#' @return a list of blocks per segment 
#' @export

endreCluster <- function(bb=bb, bwd=300, min.block=3, hclust=0.4)
        {
        pos    <- as.integer(unlist(lapply(strsplit(rownames(bb),split="\\."), "[[", 2)))
        d      <- methylDist(bb,pos, bwd = bwd)
        clRes  <- hclust(as.dist(d),method = 'single')
        cutRes <- cutree(clRes,h = hclust)
        tab    <- table(cutRes)
        temp   <- aggregate(names(cutRes), by=list("members"=as.integer(cutRes)), c)
        temp   <- temp[temp$member %in% names(which(tab >=min.block)),]$x
        if(is.matrix(temp)) {temp <- list(as.vector(temp))}
        if(length(temp) == 0){return(c())}
        return(temp)
        }

#' getBlocks
#'
#' Iterates over all segments 
#' Divides all chromsome segments into blocks based on covariation and distance
#' Limits index to blocks
#' Blocks are stored in object@blocks
#' Updated index stored in object@index
#' @param object an object of meth class 
#' @param bwd bandwidth
#' @param min.block minumum number CpGs sites per block
#' @param hclust hclust cutoff for dividing blocks
#' @param ncores number of cores to use 
#' @return a list of blocks for all segments on a chromsome
#' @export

getBlocks <- function(object, bwd=300, min.block=5, hclust=0.4, ncores=16)
        {
	#### paramaters ####
	b0    <- object@imp_mat
        achr  <- object@achr
	segs  <- object@segs
	index <- object@index	
        #### find blocks from segments ####     
        t0        <- mclapply(1:length(segs), function(i){endreCluster(bb=b0[segs[[i]],], bwd=bwd, min.block=min.block, hclust=hclust)}, mc.cores=ncores)
        names(t0) <- names(segs)
        t0        <- unlist(t0, recursive=F)
        t0        <- t0[which(as.numeric(lapply(t0, length)) > 0)]
        if (length(t0) < 1)
		{
		object@blocks <- GRanges() 
		object@index  <- GRanges() 
		object@imp_mat <- data.frame()  
		return(object) 
		}		
	### what if there are no blocks ###
        ### create block GRanges and update index###
        blocks         <- getGRangeIndex(t0, achr)
        names(blocks)  <- paste0(achr,".", 1:length(blocks))
        ov             <- findOverlaps(index, blocks)
        index$block    <- NA
        index[queryHits(ov)]$block <- names(blocks[subjectHits(ov)])
        index <- index[!is.na(index$block),]
	#### return ####
	object@blocks <- blocks 
	object@index  <- index 
	return(object) 
        }

#' filterIQR
#'
#' Filters blocks whose maximum-minimum divided by iqr is greater than cutoff
#' This removes blocks were one or two sites are responsible for variation
#' Updated index stored in object@index
#' Updated blocks stored in object@blocks 
#' @param object an object of meth class 
#' @param iqr_cutoff the iqr cutoff 
#' @return an object of meth class
#' @export

filterIQR <- function(object, iqr_cutoff=10)
        {
	if (length(object@blocks) < 1){return(object)}
	#### paramaters ####
	b0     <- object@imp_mat 
	index  <- object@index 
	blocks <- object@blocks 
	##### remove outliers ####	
        temp_imp      <- b0[names(index),]
        xx            <- split(temp_imp, index$block)
        xx            <- lapply(xx, as.matrix)
        iqr           <- lapply(xx, IQR)
        rr            <- lapply(lapply(xx, range), diff)
        fx            <- as.numeric(rr)/as.numeric(iqr[names(rr)])
        fx            <- setNames(fx, names(xx))
        fx            <- fx[is.finite(fx)]
        fx            <- fx[which(fx <= iqr_cutoff)]
	keep          <- names(fx) 
	#### check to see if we have no blocks 
	if (length(keep) < 1) 
		{
		object@blocks <- GRanges() 
		object@index <- GRanges() 
		}
	if (length(keep)>0)
		{	
		object@blocks <- sort(blocks[keep])
        	object@index  <- index[index$block %in% keep]
		}
	return(object) 
        }

#' calcBlockCoverage
#'
#' calculates block statistics for analysis
#' calculates average number of CpGs per block (block_mean)
#' calculates total number reads per block (block block_sum)
#' calculates number of CpGs per block (L) 
#' calculates average methylation level per block (PC)
#' calculates average impution per block (imp) 
#' @param object an object of meth class 
#' @param min.block minimum number of CpG sites per block 
#' @param ncores number of cores to use 
#' @return an object of meth class
#' @export

calcBlockCoverage <- function(object, min.block=5, ncores=16) 
	{
	#### wei error ###
	if (length(object@blocks) < 1){return(object)}
	#### paramaters ####
	blocks          <- object@blocks
        index           <- object@index
	b0              <- object@imp_mat 
	cov             <- object@raw_cov
	b0              <- b0[names(index),]
	cov             <- cov[names(index),]
	cov[is.na(cov)] <- 0 
	#### initialize ####
 	blocks$block_mean <- 0
	blocks$block_sum  <- 0
	blocks$L          <- 0
	blocks$PC	  <- 0 
	blocks$imp        <- 0 
	#### calculate average methylation ###
	xx <- split(b0, index$block) 
	PC <- lapply(lapply(xx, colMeans), mean)
	blocks[names(PC)]$PC <- as.numeric(PC)	
	#### calculate number of cpgs  ####
	xx <- split(cov, index$block) 
	L  <- lapply(xx, nrow) 
	blocks[names(L)]$L <- as.numeric(L) 
	#### calculate total reads per block ####
	block_sum <- lapply(xx, sum)
	blocks[names(block_sum)]$block_sum <- as.numeric(block_sum)  
	#### calculate average reads per block ####
	block_mean <- lapply(lapply(xx, rowSums), mean) 
	blocks[names(block_mean)]$block_mean <- as.numeric(block_mean) 
	#### calculate impution ####
	imp <- lapply(xx, function(xx){sum(xx==0)/(nrow(xx)*ncol(xx))})	
	blocks[names(imp)]$imp <- as.numeric(imp) 
	### verify and update####
	blocks <- blocks[which(blocks$L >= min.block)]
        index  <- index[index$block %in% names(blocks)]
        ### sort and update opject 
	index <- sort(index) 
	object@imp_mat <- b0[names(index),]
	object@raw_cov <- cov[names(index),]
	object@blocks  <- blocks
        object@index   <- index
        return(object)
	}

#' getReturnData
#'
#' Returns list of m2=imputed matrix, blocks=blocks, index=index, raw_cov=raw_cov and segments=segs 
#' @param object an object of meth class 
#' @return list of m2=imputed matrix, blocks=blocks, index=index, raw_cov=raw_cov and segments=segs  
#' @export

getReturnData <- function(object)
        {
	### what if there are no blocks 	
        b0               <- object@imp_mat
	blocks           <- object@blocks 
	index            <- object@index
	raw_cov          <- object@raw_cov
	segs             <- object@segs
	cat(paste("found", length(blocks), "blocks\n"))  
	return(list(m2=b0, blocks=blocks, index=index, raw_cov=raw_cov, segs=segs)) 	
        }



