library(GenomicRanges)

load(rws,"~/path/results/rws_ENCODE.Rdat")
#load(rws_iid,"~/path/results/rws_ENCODE.Rdat")
load(rws_sd,"~/path/results/rws_ENCODE.Rdat")
#load(rws_sd_iid,"~/path/results/rws_ENCODE.Rdat")

regSizes <- c(24897L, 24221L, 19825L, 19023L, 18150L, 17077L, 15936L, 
              14510L, 13836L, 13381L, 13510L, 13329L, 11438L, 10691L, 10201L, 
              9025L, 8327L, 8029L, 5863L, 6436L, 4672L, 5083L, 15605L)

names(regSizes) <- paste("chr",c(1:22,"X"),sep="")

regionGR <- list()
for(k in paste("chr",c(1:22,"X"),sep=""))
{
  minpos <- 0
  reg_size <- 1e04
  maxpos <- regSizes[k]*reg_size
  nregions <- ceiling((maxpos-minpos) / reg_size)
  regstarts <- seq(minpos,(nregions-1)*reg_size+minpos,reg_size)+1
  regends <-  seq(minpos+reg_size,nregions*reg_size+minpos,reg_size)
  regionGR[[k]] <- GRanges(seqnames=Rle(rep(k,nregions)),ranges=IRanges(regstarts,regends))
  regionGR[[k]]$id <- paste(rep(k,nregions),":",regstarts,"..",regends,sep="")
}
regionGR <- GRangesList(regionGR)
regGR <- unlist(regionGR)

names.list <- colnames(rws)

rwdat <- list()
for(p in 1:length(names.list))
{
  print(p)
  differ0 <- rws[,p]
  differ1 <- c(0,diff(rws[,p]))
  sd0 <- rws_sd[,p]
  
  rw_stab <- apply(rws,1,sd)
  
  rwdat[[p]] <- data.frame("rw"=differ0,"rw_diff"=differ1,"rw_sd"=sd0,
                           "rw_stab"=rw_stab,"rw_sd_stab"=rw_sd_stab)
}

source("getRankedBounds.R")

bl <- list()
for(q in 1:length(names.list))
{
  p <- grep(paste(names.list,collapse="|"),colnames(rw_mat))[q]
  print(paste(colnames(rws)[p]))
  pred.b <- getRankedBounds(mydata=rwdat[[p]],total.aim = 5000,wq=1)
  bl[[q]] <- pred.b
}