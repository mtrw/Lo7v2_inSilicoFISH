setwd("/data/gpfs/projects/punim1869/users/mrabanuswall/workspace/rye_fish")
source("scripts/0_setup.R")

# Colours!
L1col <- "#77050577" # Colour for library 1 probes ...
L2col <- "#05058877"


# Load data (if absent, run 1_run_blast.R and 2_count_gBins.R)
# WARNING, CURRENT SAVED VERSION INCOMPLETE!!! DELETE THIS TEXT ONLY WHEN IT'S REALLY REALLY DONE
gBins <- readRDS(gBinsRdsFname)
gBins[,.N,by=.(genome)]

# readFai(refFaiFnames[1])[,sum(seqLength)/1e9]
# readFai(refFaiFnames[2])[,sum(seqLength)/1e9]
# readFai(refFaiFnames[3])[,sum(seqLength)/1e9]


# Sanity check on data # (b)in (s)ummary
bs <- gBins[ , .(totBS=sum(sumBitscore)) , by=.(genome,seqId,lib) ] %>% dcast(formula=genome+seqId~lib,value.var = "totBS")
bs

# Establish table of what genomes / chrs want plotting
plotLooper <- bs[seqId %!in% c("0R","GWHASIY00000008","chrUn"),.(g=genome,sid=seqId)] # In case you want to set order here? Keep it in order by genome or things will be ugly

# Layout
dev.off()
par(mfcol=c(plotLooper[,.N,by=.(g)][,max(N)],plotLooper[,nu(g)]),mar=c(0,1,0,1),oma=c(3,1,1,1))

# Plot!
for(pl_i in 1:nrow(plotLooper)){
  #pl_i <- 8
  p <- gBins[genome==plotLooper[pl_i]$g & seqId==plotLooper[pl_i]$sid]
  xaxtArg <- if(pl_i%%7==0){NULL}else{"n"}
  null_plot(range(p$start),range(0,p$sumBitscore*1.4),xaxt=xaxtArg,yaxt="n") # consider scaling to whole dataset ...
  lines(p[lib=="L1"]$start,p[lib=="L1"]$sumBitscore,col=L1col,lwd=2)
  lines(p[lib=="L2"]$start,p[lib=="L2"]$sumBitscore,col=L2col,lwd=2)
  mtext(plotLooper[pl_i]$sid,side=2,cex=.5)
  if(pl_i%%7==1){ mtext(plotLooper[pl_i]$g,side=3,cex=.5,padj=2) }
}






























