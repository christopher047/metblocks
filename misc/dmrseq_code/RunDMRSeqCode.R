library(data.table) 
library(GenomicRanges) 
library(bsseq)
library(dmrseq)
library(metblocks) 

set.seed(1024) 

############## read data ####################
path         <- file.path(system.file("extdata", package="metblocks"), "")
achr         <- "chr18"
tot_cov      <- read.table(paste0(path, achr, "_tot_cov.tab"), header=T, row.names=1, sep="\t") 
tot_meth     <- read.table(paste0(path, achr, "_tot_meth.tab"), header=T, row.names=1,  sep="\t") 

############ no NAs allowed set to zero #####
tot_cov[is.na(tot_cov)] <- 0 	
tot_meth[is.na(tot_meth)] <- 0 

########## create index for GRanges ####
chr          <- gsub("\\..*$","", rownames(tot_cov))
pos          <- gsub("^.*\\.","", rownames(tot_cov))
index        <- GRanges(data.frame("chr"=chr, "start"=pos, "end"=pos))
names(index) <- rownames(tot_cov) 
index        <- sort(index) 

############ make sure they are in the same order ##########
tot_cov  <- tot_cov[names(index),]
tot_meth <- tot_meth[names(index),]
tot_pct  <- tot_meth/tot_cov 

########### pheno ######
cond <- unlist(lapply(strsplit(colnames(tot_meth), "_"), "[[",1))
ids  <- unlist(lapply(strsplit(colnames(tot_meth), "_"), "[[",2))
pheno <- data.frame("id"=ids, "cond"=cond) 
rownames(pheno) <- colnames(tot_meth) 

########## run dmr seq with defaults #########
bs             <- BSseq(gr=index, M=as.matrix(tot_meth), Cov=as.matrix(tot_cov)) 
pData(bs)$cond <- pheno$cond
bs             <- sort(bs)
regions        <- dmrseq(bs=bs, testCovariate="cond")

############## results  ###############
rawDiff      <- meanDiff(bs, dmrs=regions, testCovariate="cond")
regions$rawDiff <- rawDiff
regions$Comp    <- "UC-N"
regions$Class   <- ""
hyper <- which(regions$rawDiff > 0) 
hypo  <- which(regions$rawDiff < 0) 
if (length(hyper) > 0){regions[hyper]$Class <- "hyper"}
if (length(hypo)  > 0){regions[hypo]$Class  <- "hypo"}
write.table(data.frame(regions), "./results/chr18_regions.tab", sep="\t")



