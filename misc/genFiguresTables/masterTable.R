library(Gviz) 
library(GenomicRanges) 
library(ChIPpeakAnno) 

### et qvalue for metiline and dmrseq
qvalue <- 0.1

### get block processed data ###
mat    <- read.table("../results/raw_mat_trim.tab", header=T, row.names=1, sep="\t")
index  <- GRanges(read.table("../results/raw_index_trim.tab", header=T, row.names=1, sep="\t")) 
blocks <- GRanges(read.table("../results/raw_blocks_trim.tab", header=T, row.names=1, sep="\t")) 

### get metilene data ####
mets      <- read.table("../metilene_code/formatted/chr18_formatted.tab", header=T, row.names=1, sep="\t")  
mets      <- GRanges(mets) 
sig_mets  <- mets[mets$qvalue <= qvalue]

### get dmrseq data ###
dmr_file  <- read.table("../dmrseq_code/results/chr18_regions.tab", header=T, row.names=1, sep="\t")
dmrs      <- GRanges(dmr_file)
dmrs      <- sort(dmrs)
sig_dmrs  <- dmrs[dmrs$qval <= qvalue]

### get total sites in all chromosomes ### 

cov_file <- read.table(system.file("extdata/chr18_tot_cov.tab", package="metblocks"), header=T, row.names=1, sep="\t")
met_file <- read.table(system.file("extdata/chr18_tot_meth.tab", package="metblocks"), header=T, row.names=1, sep="\t")
tot_sites <- sum(!is.na(cov_file)) 

### table basic sites 1###
r1             <- c(length(blocks), length(sig_dmrs), length(sig_mets)) 
### percentage of total sites in result
r2             <- c(sum(blocks$L)/tot_sites, sum(sig_dmrs$L)/tot_sites, sum(sig_mets$num_CpGs)/tot_sites)*100   
### average block size 
r3             <- c( mean(width(blocks)), mean(width(sig_dmrs)), mean(width(sig_mets)))
### average number of CpGs in block 
r4             <- c(mean(blocks$L), mean(sig_dmrs$L), mean(sig_mets$num_CpGs))
tab1           <- do.call("rbind", list(r1, r2, r3, r4))
colnames(tab1) <- c("blocks", "dmrseq", "metilene") 
rownames(tab1) <- c("number", "pct of total", "mean width", "mean number of CpGs")
tab1           <- round(tab1, 2) 
write.table(tab1, "../tables/table1.tab", sep="\t") 

#### show examples of overlaps  #####

### load meth pct from raw_files###

chromosome <- gsub("\\..*", "", rownames(cov_file))
end        <- gsub(".*\\.", "", rownames(cov_file))
raw_pct    <- GRanges(data.frame("chromosome"=chromosome, "start"=end, "end"=end, met_file/cov_file))
ov         <- findOverlaps(blocks, sig_dmrs)
temp       <- aggregate(queryHits(ov), list(subjectHits(ov)), c)
temp       <- temp[order(unlist(lapply(temp$x, length)), decreasing=T),]
all_exons  <- GRanges(read.table(system.file("extdata/exons_chr18.tab", package="metblocks"), sep="\t"))  

#### start loop here ###
for(ii in 1:5) 
	{
	#### data track blocks ####
	tdat       <- index[index$block %in% names(blocks[temp[ii,2][[1]]])] 
	nn         <- grep("^NN", colnames(mat)) 
	uc         <- grep("^UC", colnames(mat)) 
	tdat$nn    <- as.numeric(rowMeans(mat[names(tdat),nn]))  
	tdat$uc    <- as.numeric(rowMeans(mat[names(tdat),uc]))  
	tdat       <- tdat[,c("nn", "uc")]
	tbs        <- blocks[temp[ii,2][[1]]]

	#### data track full region  ####
	tpct       <- GRanges(paste0(unique(as.character(seqnames(tdat))) , ":", min(start(tdat)), "-", max(end(tdat))))
	tpct       <- raw_pct[queryHits(findOverlaps(raw_pct, tpct))]  
	nn         <- grep("^NN",  colnames(mcols(tpct))) 
	uc         <- grep("^UC",  colnames(mcols(tpct))) 
	tpct$nn    <- as.numeric(rowMeans(data.frame(mcols(tpct)[,nn]), na.rm=T))
	tpct$uc    <- as.numeric(rowMeans(data.frame(mcols(tpct)[,uc]), na.rm=T))
	tpct       <- tpct[, c("nn", "uc")]

	#### make sure the positions align with start end 
	raw_track        <- DataTrack(tpct, name="raw %") 
	meth_track       <- DataTrack(tdat, name="blocks imp %")
	axis_track       <- GenomeAxisTrack()
	block_track      <- GeneRegionTrack(blocks[temp[ii,2][[1]]], name="metblocks") 
	dmrseq_track     <- GeneRegionTrack(sig_dmrs[temp[ii,1]], name="dmrseq") 
	met_track        <- GeneRegionTrack(sig_mets[unique(queryHits(findOverlaps(sig_mets, tbs)))], name="metilene") 
	exon_track       <- GeneRegionTrack(all_exons[unique(queryHits(findOverlaps(all_exons, tbs)))], name="exon")


	displayPars(raw_track) <-  list("groups"=c("nn", "uc"), "col"=c("green", "red"), "col.frame"="lightgrey", "col.axis"="black") 
	displayPars(raw_track)  <- list("box.legend"=TRUE, "cex.legend"=1) 
	displayPars(meth_track) <- list("groups"=c("nn", "uc"), "col"=c("green", "red"), "col.frame"="lightgrey", "col.axis"="black") 
	displayPars(meth_track) <- list("legend"=FALSE) 

	all_tracks <- list(raw_track, meth_track, axis_track, block_track, dmrseq_track, met_track, exon_track)
	names(all_tracks) <- c("raw_track","meth_track","axis_track","block_track","dmrseq_track","met_track","exon_track")
	scheme     <- list("rotation.title"=0, "fontcolor.title"="black", "cex.title"=1,  "background.title" = "white", "frame"=TRUE) 

	for(i in 1:length(all_tracks)){displayPars(all_tracks[[i]]) <- scheme}
	displayPars(all_tracks$raw_track) <- list("rotation.title"=90) 
	displayPars(all_tracks$meth_track) <- list("rotation.title"=90) 
	displayPars(all_tracks$axis_track) <- list("frame"=FALSE) 
	displayPars(all_tracks$exon_track) <- list("stacking"="squish") 

	### remove empty tracks ###
	keep <- c(names(all_tracks)[1:3], names(which(lapply(all_tracks[4:length(all_tracks)], length) > 0)))
	all_tracks <- all_tracks[keep]

	pdf(paste0("../figures/test_", sprintf("%02d", ii), ".pdf"), width=10, height=10)  
	plotTracks(all_tracks, from=min(start(tdat))-5, to=max(end(tdat))+5, title.width=2)  
	dev.off() 
	#### end loop
	cat(ii, sep="\n")  
	}

### density plot of regions ####
pdf("../figures/region_density.pdf", width=6, height=6) 
plot(density(width(blocks)), xlim=c(0,800), xlab="Width", main="", cex=1.5)  
lines(density(width(sig_dmrs)), col="red") 
lines(density(width(sig_mets)), col="blue") 
legend("topright", legend=c("metblocks", "metiline", "dmrseq"), col=c("black", "blue", "red"), lty=1)
dev.off() 

### show overlaps in a venn ###
pdf("../figures/venn_overlap_regions.pdf", width=6, height=6)
venn <- makeVennDiagram(list(blocks, sig_mets, sig_dmrs), c("metblocks", "metilene", "dmrseq"), by="region", col=c("black", "blue", "red"), maxgap=10, cex=1.5)
dev.off()

###pca plot of seperation###
m2 <- read.table("../results/raw_mat_trim.tab", header=T, row.names=1, sep="\t") 
cols <- rep("green", ncol(m2)) 
cols[grep("UC", colnames(m2))] <- "red"
dd <- prcomp(t(m2)) 
pc1 <- paste0(round(summary(dd)$importance[2,1]*100,1), "%") 
pc2 <- paste0(round(summary(dd)$importance[2,2]*100,2), "%") 
xlab <- paste("PC1", pc1) 
ylab <- paste("PC2", pc2) 

pdf("../figures/pca_plot_chr18.pdf", width=6, height=6) 
plot(dd$x, col=cols, pch=19, cex=2, main="", xlab=xlab, ylab=ylab)
dev.off() 

### move to article ###
file.rename("../tables/table1.tab", "../article/table1.tab") 
file.rename("../figures/venn_overlap_regions.pdf", "../article/venn_overlap_regions.pdf") 
file.rename("../figures/test_05.pdf", "../article/example_overlap.pdf") 
file.rename("../figures/region_density.pdf", "../article/region_density.pdf") 
file.rename("../figures/pca_plot_chr18.pdf", "../article/pca_plot_chr18.pdf") 









