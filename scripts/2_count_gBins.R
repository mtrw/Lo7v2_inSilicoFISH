setwd("/data/gpfs/projects/punim1869/users/mrabanuswall/workspace/rye_fish")
source("scripts/0_setup.R")

# Read blast files and aggregate THIS WILL ITERATIVELY UPDATE `gBins`. Do it after a fresh run of `0_setup.R`.
# Requires blast output-->if absent run 1_run_blast.R
{
  l_ply(unique(probeSeqs$blGrp),function(bg){
    l_ply(1:length(refFastaFnames),function(rfn_i){
      #dev bg <- 1; rfn_i <- 3
      blFname <- paste0(blastDir,"/blastProbes_gp:",bg,"_genome:",refGenomeNames[rfn_i],".rds")
      if(!fileExistsNonZeroSize(blFname)){
        stop("ERROR: Non-existent file ",blFname,"\nAfter stopping like this, be sure to reinitialise `gBins` before restarting!")
      } else {
        #browser()
        bl <- readRDS(blFname)
        olBl <- bl[,.( # for overlapping
          seqId=sSeqId,
          start=tmp<-round(pmean2(sStart,sEnd)),
          end=tmp,
          lib=qSeqId,
          genome,
          bitscore
        ),] %>% setkey(lib,genome,seqId,start,end)
        addMe <- foverlaps(gBins,olBl,nomatch=0L)[,.(countAdd=.N,bitscoreAdd=sum(bitscore)),by=.(lib,genome,seqId,i.start,i.end)] # find bins that overlap, for each bin, total up the counts and bitscores
        gBins[addMe,c("count","sumBitscore"):=.( # ... and add these to the table of genomic bins
          count+countAdd,
          sumBitscore+bitscoreAdd
        )]
        ce(blFname," added to bin scores ...")
      }
      NULL
    }) %>% invisible
    NULL
  }) %>% invisible
  
  
  
  saveRDS(gBins,gBinsRdsFname)
}



