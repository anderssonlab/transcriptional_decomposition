getRankedBounds <- function(mydata=rwdat[[1]],total.aim=5000,filter=FALSE,filter_size=10,wq=0.5,myreg=regGR)
{
  chr_aim <- list()
  
  seq.lengths <- tapply(as.vector(seqnames(myreg[which(mydata[,"rw"]>0)])),as.vector(seqnames(myreg[which(mydata[,"rw"]>0)])),length)
  bounds.per.chr <- round(total.aim*(seq.lengths/sum(seq.lengths)))
  
  for(k in paste("chr",c(1:22,"X"),sep=""))
  {
    mybins <- (which(as.vector(seqnames(myreg))==k))
    myw <- mydata[mybins,"rw"]
    mysd <- mydata[mybins,"rw_sd_stab"]
    aim <- bounds.per.chr[k]
    mydiff <- (mydata[mybins,"rw_diff"])
    myd <- abs(mydiff)/mysd
    myd <- (myd * (wq*myw))
    
    crosses <- extrema(myd)$maxindex[,1]
    summary(abs(myd[crosses]))
    acrosses <- sort(crosses[order(abs(myd[crosses]),decreasing=T)[1:aim]])
    summary(abs(myd[acrosses]))
    summary(diff(acrosses))
    
    if(filter)
    {
      while( length(which(diff(acrosses)<filter_size))>0 )
      {
        use <- which(diff(acrosses)<filter_size)[sample(1:length(which(diff(acrosses)<filter_size)),1)]
        mywdc <- abs(myd)[acrosses[(use):(use+1)]]
        acrosses[use:(use+1)] <- acrosses[(use):(use+1)][which(mywdc==max(mywdc))]
        acrosses <- unique(acrosses)
      }
    }
    
    chr_aim[[k]] <- mybins[acrosses]
    print(paste(k,length(chr_aim[[k]])))
  }
  
  ab <- unlist(chr_aim)
  
  myblist <- list()
  myblist[["up"]] <- rep(0,length(myreg))
  myblist[["up"]][ab[mydata[,"rw_diff"][ab]>0]] <- 1
  myblist[["down"]] <- rep(0,length(myreg))
  myblist[["down"]][ab[mydata[,"rw_diff"][ab]<0]] <- 1
  lapply(myblist,sum)
  
  ab <- sort(c(which(myblist[["up"]]==1)+1,which(myblist[["down"]]==1)))
  
  return(ab)
}