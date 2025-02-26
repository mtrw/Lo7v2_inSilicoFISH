setwd("/data/gpfs/projects/punim1869/users/mrabanuswall/workspace/rye_fish")
source("scripts/0_setup.R")

# Want to start from scratch? Run me
#unlink(paste0(blastDir,"/*"))


# Run blast
require(parallel)
mcldtply(mc.cores=5,unique(probeSeqs$blGrp),function(bg){
  mcldtply(mc.cores=3,1:length(refFastaFnames),function(rfn_i){
    #dev bg <- 1; rfn_i <- 3
    blFname <- paste0(blastDir,"/blastProbes_gp:",bg,"_genome:",refGenomeNames[rfn_i],".rds")
    if(!fileExistsNonZeroSize(blFname)){
      saveRDS(blast(querySeqDT=probeSeqs[blGrp==bg],subjectFname=refFastaFnames[rfn_i],moreBlastArgs=blastExtraArgs)[,genome:=refGenomeNames[rfn_i]],file=blFname)
    }
    NULL
  }) %>% invisible
  NULL
}) %>% invisible



