##required inputs

## CAGE library IDs for  sample replicates
#p=1 refers to the first cell type of the random walk matrix (i.e. GM12878)
#enhancers is the matrix of enhancer expression with corresponding enhancerGR (locations of corresponding enhancers)
#groupGR is the dpi cluster locations
#sample.names is the names of the samples
#use_CTCF - include the bin's CTCF information?
#myctcf - vector of length equal to the number of bins with a number representing the number of significantly bound CTCF sites within that bin
#myiid - iid vector
#rws_sd_iid=iid_sd_mat
#rws_sd - random walk standard deviation
#rws_sd_iid_stab - independent component stability
#rws_sd_stab - random walk stability
#myrws - random walk
#boundSet - set of boundaries, a vector of bin numbers which are designated boundary bins
#myDIST - maximum distance between pairs to consider
#gmMat - raw expression matrix of called promoters (dpi clusters)
#includeGMcounts - whether or not to include the interaction count in the dataset outputed
#mychrs - which chromosomes to calculate interaction dataset over
#peaksGR - locations of dpi clusters
#myregGR - locations of 10kb bins
#toNORM - whether or not to normalise outputed datasets
#EN_ONLY - only consider enhancer-promoters interactions (not promoter-promoter)

getTDAT <- function(p=1,p_libs=c("CNhs12331", "CNhs12332", "CNhs12333"),enhancers=enhancers,groupGR=ENCdpiGR,sample.names=samp.names,use_CTCF=FALSE,myctcf=CTCF_vec,myiid=iid_mat,rws_sd_iid=iid_sd_mat,rws_sd=rw_sd_mat,rws_sd_iid_stab=rws_sd_iid_stab,rws_sd_stab=rws_sd_stab,myrws=rw_mat,boundSet=bl[[1]],myDIST=2e06,gmMat=gmMat,includeGMcounts=FALSE,mychrs=paste("chr",c(1:22,"X"),sep=""),peaksGR=peaksGR,myregGR=regGR,toNORM=FALSE,EN_ONLY=TRUE)
{
  print(paste("collecting enhancer and promoter information..."))
  
  en_sup <- rowSums(enhancers_bin[rowSums(enhancers)>0,])
  
  enhancers2 <- enhancers[,p_libs]
  
  orth_h <- row.names(enhancers)
  chr_h <- gsub("^(chr.+):.+$","\\1",orth_h)
  st_h <- gsub("^.+:(.+)-.+$","\\1",orth_h)
  en_h <- gsub("^.+:.+-(.+)$","\\1",orth_h)
  len_h <- as.numeric(en_h)-as.numeric(st_h)+1
  enhancerGR <- GRanges(seqnames=Rle(chr_h),ranges=IRanges(as.numeric(st_h),as.numeric(en_h)),name=orth_h,length=len_h)
  
  peaksGR <- GRanges(Rle(as.vector(peaks$chr)),IRanges(as.vector(peaks$start),as.vector(peaks$end)),name=peaks$name,gene=peaks$gene_name)
  matLibs <- gmMat[,as.vector(p_libs)]
  matLibs <- matLibs[which(apply(matLibs,1,sum)>0),]
  mylds <- colSums(matLibs)
  matLibs <- matLibs[which(apply(matLibs,1,function(x) length(x[x>0]))>1),]
  
  for(erow in 1:ncol(enhancers2))
  {
    enhancers2[,erow] <- (enhancers2[,erow]/1e04)/(mylds[erow]/1e04/1e06)
  }
  eRNA <- rowMeans(enhancers2)
  print(sum(eRNA))
  
  
  print(paste("adding promoter and enhancer information..."))
  
  #number of peaks in GM
  nopeaks <- countOverlaps(myregGR,peaksGR[peaksGR$name %in% row.names(matLibs)])
  
  #number of enhancers in GM
  noens <- countOverlaps(myregGR,enhancerGR[eRNA>0])
  ensup <- countOverlaps(myregGR,enhancerGR)

  #mean number of eRNAs
  enRNA <- rep(0,length(myregGR))
  ensup_ol <- findOverlaps(myregGR,enhancerGR)
  enRNA[unique(queryHits(ensup_ol))] <- tapply(eRNA[subjectHits(ensup_ol)],queryHits(ensup_ol),mean)
  
  #get predicted boundaries from the random walk
  boundList <- list()
  
  boundList[["moves"]] <- rep(0,nrow(myrws))
  boundList[["moves"]][boundSet] <- 1
  
  mye1 <- e1vec
  
  mybvec <- boundList[["moves"]]
  mybvec <- cumsum(mybvec) + 1
  
  mybounds2 <- extrema(e1vec)$cross[,1]
  mybvec2 <- rep(0,length(myregGR))
  mybvec2[mybounds2] <- 1
  mybvec2 <- cumsum(mybvec2) +1 
  
  mytrain_list_data <- list()
  for(k in mychrs)
  {
    print(paste("processing...",k))
    
    if(includeGMcounts)
    {
      print(paste("opening hic data...",k))
      
      int <- ints_list[[k]]
      starts <- seq(0,(nrow(int)-1)*10000,by=10000)
      ends <- seq(10000,(nrow(int))*10000,by=10000)
      chr <- rep(k,nrow(int))
      hicGR <- GRanges(seqnames=Rle(chr),ranges=IRanges(as.numeric(starts)+1,as.numeric(ends)),id=c(1:length(chr)))
      
      myintuse <- which(countOverlaps(hicGR,myregGR)>0)
      mytouse <- which(countOverlaps(myregGR,hicGR)>0)
      int <- int[myintuse,myintuse]
    }else{
      mytouse <- which(as.vector(seqnames(myregGR))==k)
    }
    
    gm_plus <- rowSums(mymat_plus[,sample.names])
    gm_minus <- rowSums(mymat_minus[,sample.names])
    
    mydir <- (gm_plus-gm_minus)/(gm_plus+gm_minus+1)
    
    rwdat_use <- cbind(myrws[mytouse,p],rws_sd[mytouse,p],rws_sd_stab[mytouse],mydir[mytouse])
    colnames(rwdat_use) <- c("rw","rw_sd","rw_sd_stab","dir")
    
    #CTCF
    if(use_CTCF==TRUE)
    {
      rwdat_use <- cbind(rwdat_use,myctcf[mytouse]) 
      colnames(rwdat_use) <- c("rw","rw_sd","rw_sd_stab","dir","CTCF")
    }
    
    iiddat_use <- cbind((myiid[mytouse,p]),(rws_sd_iid[mytouse,p]),rws_sd_iid_stab[mytouse])
    colnames(iiddat_use) <- c("iid","iid_sd","iid_sd_stab")
    
    mytouse_chr <- which(countOverlaps(myregGR[as.vector(seqnames(myregGR))==k],myregGR[mytouse])>0)
    
    print(paste("getting cross correlations..."))
    
    cc_mat <- cross_cor[[k]][as.vector(mytouse_chr),as.vector(mytouse_chr)]
    cc_mat_iid <- cross_cor_iid[[k]][as.vector(mytouse_chr),as.vector(mytouse_chr)]
    
    npeaks <- nopeaks[mytouse]
    nens <- noens[mytouse]
    en_support <- ensup[mytouse]
    eRNAsum <- enRNA[mytouse]
    e1v <- mye1[mytouse]
    
    dat_use <- cbind(rwdat_use,iiddat_use,npeaks,nens,en_support,eRNAsum,e1v)
    
    print(paste("calculating pairs..."))
    
    
    mypairs <- as.data.frame(expand.grid(1:length(mytouse),1:length(mytouse)))
    colnames(mypairs) <- c("bait","targ")
    mypairs$bait_id <- mytouse[as.vector(mypairs[,1])]
    mypairs$targ_id <- mytouse[as.vector(mypairs[,2])]
    mypairs$dist <- abs( mypairs$bait_id -  mypairs$targ_id )
    
    mypairs <- mypairs[mypairs$dist <= (myDIST/1e04),]
    
    mygenebaits <- which(countOverlaps(regGR[mypairs$bait_id],mybaits.GR)>0)
    mypairs <- mypairs[mygenebaits,]

    if(EN_ONLY)
    {
      print(paste("reducing targets to enhancers only..."))
      
      myentargs <- which(countOverlaps(myregGR[mypairs$targ_id],enhancerGR)>0)
      mypairs <- mypairs[myentargs,]
    }else{
      print(paste("reducing targets..."))
      
      mygenetargs <- which(countOverlaps(myregGR[mypairs$targ_id],groupGR)>0)
      length(mygenetargs)
      mypairs <- mypairs[mygenetargs,] 
      print(dim(mypairs))
      
      mypairs <- mypairs[-which(mypairs$targ_id==mypairs$bait_id),] 
      
      torm <- (paste(mypairs$targ,mypairs$bait) %in% paste(mypairs$bait,mypairs$targ))
      mypairs <- mypairs[-which(torm),] 
    }
    
    mypairs$nbounds <- apply(mypairs[,3:4],1,function(x) abs(mybvec[x[1]] - mybvec[x[2]]) )
    mypairs$nbounds[mypairs$nbounds>5] <- 5 
    
    mypairs$eigcross <- apply(mypairs[,3:4],1,function(x) abs(mybvec2[x[1]] - mybvec2[x[2]]) )
    mypairs$eigcross[mypairs$eigcross>5] <- 5
    
    mypairs$en_bet <- apply(mypairs[,3:4],1,function(x) sum(enRNA[x[1]:x[2]]) - sum(enRNA[c(x[1],x[2])]) )
    
    mypairs$bait_in_group <- countOverlaps(myregGR[mypairs$bait_id],groupGR)
    
    d1 <- dat_use[mypairs[,1],]
    d2 <- dat_use[mypairs[,2],]
    colnames(d1) <- paste("bait_",colnames(d1),sep="")
    colnames(d2) <- paste("targ_",colnames(d2),sep="")
    
    myinfo <- as.data.frame(cbind(mypairs,d1,d2)) 
    
    myinfo$rwdiff <- abs(myinfo$bait_rw - myinfo$targ_rw)
    myinfo$iiddiff <- abs(myinfo$bait_iid - myinfo$targ_iid)
    
    if(toNORM)
    {
      myinfo[,-c(1:6)] <- apply(myinfo[,-c(1:6)],2,function(x) scale(x)[,1])
    }
    
    myinfo$cross_rw <- as.vector(cc_mat[as.matrix(mypairs[,1:2])])
    myinfo$cross_iid <- as.vector(cc_mat_iid[as.matrix(mypairs[,1:2])])
    myinfo$chr <- rep(k,nrow(myinfo)) 
    
    if(includeGMcounts)
    {
      myinfo$count <- as.vector(int[as.matrix(mypairs[,1:2])])
    }
    
    mytrain_list_data[[k]] <- myinfo
    print(paste(dim(myinfo)))
  }
  return(mytrain_list_data)
}

