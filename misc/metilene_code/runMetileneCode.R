
### metilene must be in path ###

#### os command to run metilene####
cmd  <- file.path(getwd(), "run_metiline.cmd") 
system(cmd) 

### write results ###
in_files <- list.files(path="./results", full.name=T) 
out_files <- paste0("./formatted/", gsub("_met_result", "_formatted", basename(in_files)))
for(i in 1:length(in_files)) 
	{
	fs <- file.info(in_files[i])$size[1]
	if(fs < 1){next}
	m2 <- read.table(in_files[i], sep="\t") 
	colnames(m2) <-c("chr","start","end","qvalue", "difference", "num_CpGs", "p_MWU", "p_2D KS","mean_NN", "mean_UC") 
	write.table(m2, out_files[i], quote=F, sep="\t") 
	}


