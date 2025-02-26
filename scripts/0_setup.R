# devtools::install_github("http://github.com/mtrw/BioDT")
# To ensure it installs, manually install msa, biostrings, data.table, plyr, magrittr
# Depends on `system()` access to e.g. blast, bedtools, awk, ... standard bioinfo linux package
library(BioDT)
source("scripts/0_functions.R")


# Filenames
refFastaFnames <- c(
  "../../../../shared_data/refs/Ensembl_Secale_cereale.Rye_Lo7_2018_v1p1p1/assembly/Secale_cereale.Rye_Lo7_2018_v1p1p1.dna.toplevel_pseudomolecules.fa",
  "../../../../shared_data/refs/Secale_cereale_Weining_Li2021/assembly/Secale_cereale_Weining_Li2021_pseudomolecules.fasta",
  "../../../../shared_data/refs/Secale_cereale_Lo7_v3_2024/assembly/Secale_cereale_Lo7_v3_2024_assembly.fasta"
)
refFaiFnames <- paste0(refFastaFnames,".fai")
refGenomeNames <- c("Lo7_2021","Weining_2021","Lo7_2024") # must be same length
readFai(refFaiFnames)[,sum(seqLength)/1e6,by=.(faiFname)]


# Kill chrs (chrUn etc)
killChrames <- c("0R","GWHASIY00000008","chrUn")

probeSeqsFname <- "data/probeseqs_lib1_lib2.fasta"
file.size(probeSeqsFname)/1e6 # not too big, cool to keep in memory

blastDir <- existentDir("data/blast")
f36dir <- existentDir("data/f36")

# Binaries
f36bin <- "~/bin/fasta36/bin/fasta36"

# (g)enomic (B)ins rds file to save final hit counts in
gBinsRdsFname <- "data/genomicBinCounts.rds"

# To save filtered blast results
blastFilteredRdsFname <- "data/blastFiltered.rds"

# Settings
  # genomic bins width and increment
binWidth <- 1e6L
binInc <- 0.5e6L
  # blast filter cutoffs for direct hit plotting
bitscoreMin <- 65


# Parameters
# data.table::setDTthreads() # to limit threads
  # blast
blastGroupSize <- 200L
blastExtraArgs <- "-task blastn-short"
  # fasta36



# Set up query (probe) seqs in DT
probeSeqs <- readFasta(probeSeqsFname,saveFnames=FALSE)[,seqId:=fifelse(grepl("lib1$",seqId),"L1","L2")]
probeSeqs[,blGrp:=rep(1:.N,each=blastGroupSize,l=.N)] %>% setkey(blGrp)
#probeSeqs[,blFname:=paste0(blastDir,"/grp",blGrp,"_bgs",blastGroupSize,".rds")]
probeSeqs[,nu(blGrp)] # Cool not too many


# Set up genomic bins to accumulate hit counts
t1 <- ldtply(seq_along(refFaiFnames),function(ffn_i){
  readFai(refFaiFnames[ffn_i])[,.(
    seqId,
    start=tmp<-seq(1L,seqLength-binWidth,by=binInc),
    end=tmp+(binWidth-1L),
    count=0L,
    sumBitscore=0.0
  ),by=.I][,genome:=refGenomeNames[ffn_i]]
})

# (g)enomic (B)ins
gBins <- rbindlist(list(copy(t1)[,lib:="L1"],copy(t1)[,lib:="L2"]))  %>% setkey(lib,genome,seqId,start,end) # THIS IS SHIT CODE, DON'T WRITE CODE LIKE THIS. Use expand.grid

gBins[]

rm(t1)


warning("IMPORTANT: At present this script depends upon the chromosome names in each genome's fasta file being different to each other. If this is not the case, UPDATE THE SCRIPT, e.g. by overlapping on one extra column containing the genome info!!!")
