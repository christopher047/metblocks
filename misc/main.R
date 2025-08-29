### requirements ###
#require(impute)
#require(data.table)
#require(GenomicRanges)
#require(bsseq)
#require(dmrseq) 
#require(Gviz)
#require(rtracklayer)
#require(ChIPpeakAnnno) 

### !!metilene must be in path!! ###

set.seed(1024) 
base       <- getwd()

### run blocks ###
setwd(file.path(base, "metblock_code")) 
source("runMetblockCode.R") 
setwd(base) 

### metilene ###
setwd(file.path(base, "metilene_code")) 
source("runMetileneCode.R") 
setwd(base) 

### dmrseq  ###
setwd(file.path(base, "dmrseq_code"))
source("RunDMRSeqCode.R") 
setwd(base) 

### tables ###
setwd(file.path(base, "genFiguresTables")) 
source("masterTable.R") 
setwd(base) 


