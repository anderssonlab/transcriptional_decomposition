#Transcriptional decomoposition script used to generate results for ENCODE CAGE datasets

library(GenomicRanges)
library(reshape2)
library(INLA)

#.bed files for relevant CAGE library CTSSs
filenames <- list.files("~/", pattern="*.bed$", full.names=TRUE)
data <- lapply(filenames, function(x) read.table(x,header=FALSE))

datnames <- gsub(".*/(.*).hg38.nobarcode.ctss.bed","\\1",filenames)
myld <- unlist(lapply(data,function(x) sum(x$V5)))
names(myld) <- datnames

#set up genomic range for each library
dataGR <- list()
for(i in 1:length(data))
{
  st_d <- as.numeric(data[[i]]$V2)
  en_d <- as.numeric(data[[i]]$V3)
  chr_d <- as.vector(data[[i]]$V1)
  tag_d <- as.vector(data[[i]]$V5)
  str_d <- as.vector(data[[i]]$V6)
  name_d <- as.vector(data[[i]]$V4)
  dataGR[[i]] <- GRanges(seqnames=Rle(chr_d),strand=Rle(str_d),ranges=IRanges(st_d,en_d),count=tag_d,id=name_d)
}

#set up and run transcriptional decomposition per chromosome (free parameters)
data_list <- list()
for(k in (paste("chr",c(1:22,"X"),sep="")))
{
  print(k)
  maxvec <- c()
  for(p in 1:length(data))
  {
    maxvec <- c(maxvec,max(data[[p]][data[[p]]$V1==k,]$V3))
  }
  
  minpos <- 0
  reg_size <- 1e04 #10kb binned data
  maxpos <- max(maxvec) + reg_size
  nregions <- ceiling((maxpos-minpos) / reg_size)
  regstarts <- seq(minpos,(nregions-1)*reg_size+minpos,reg_size)+1
  regends <-  seq(minpos+reg_size,nregions*reg_size+minpos,reg_size)
  regionGR <- GRanges(seqnames=Rle(rep(k,nregions)),ranges=IRanges(regstarts,regends))
  regionGR$id <- paste(rep(k,nregions),":",regstarts,"..",regends,sep="")
  
  #overlap the chunks with the genomic ranges to get the binning CAGE tags
  mymat <- matrix(nrow=length(regionGR),ncol=length(data))
  for(p in 1:length(data))
  {
    myOLs <- findOverlaps(regionGR,dataGR[[p]])
    mycounts <- sapply(split(dataGR[[p]][subjectHits(myOLs)]$count,queryHits(myOLs)),sum)
    myregvec <- rep(0,length(regionGR))
    myregvec[as.numeric(names(mycounts))] <- mycounts
    mymat[,p] <- myregvec
  }
  colnames(mymat) <- datnames
  data_list[[k]] <- mymat
}

mymat <- do.call(rbind,data_list)

mydat <- melt(mymat)

mydat$sample_stem <- gsub("^(.+)_.*$","\\1",mydat$Var2)
colnames(mydat) <- c("geneid","sampleid","value","sample_stem")

mydat$groupid <- as.numeric(factor(mydat$sample_stem))
mydat$groupid2 <- mydat$groupid 

mydat$sample <- mydat[,2]
mydat <- mydat[order(mydat$geneid),]
mydat$geneid2 <- mydat$geneid
mydat$ld <- myld[as.vector(mydat$sampleid)]/1e06
mydat$tissid <- as.numeric(factor(mydat$sample))
mydat$tissid2 <- mydat$tissid

formula <- value ~ 1 + f(geneid2,model="iid",replicate=groupid,hyper=prior.iid) + f(geneid,replicate=groupid,model="rw1",hyper=prior.rw,scale.model=T)

mytheta1  <- list(initial=-3)
mytheta2  <- list(initial= 2,fixed=FALSE,prior="gaussian",param=list(-.1,.1))

prior.rw <- list( 
  prec = list( 
    initial = log(0.000001), 
    fixed = FALSE)) 

prior.iid <- list( 
  prec = list( 
    initial = log(0.000001), 
    fixed = FALSE)) 

#two pass
result_1 = inla(formula,family="zeroinflated.nbinomial0",verbose=TRUE, data = mydat, 
                control.family = list(hyper = list(theta=mytheta1)),E=mydat$ld,
                control.fixed= list(prec.intercept = 1),
                control.inla = list(int.strategy="eb",strategy="gaussian",diagonal=1),
                working.directory = paste("~/path/scratch/",k,sep=""),
                num.threads = 35)


result_1 = inla(formula,family="zeroinflated.nbinomial0",verbose=TRUE, data = mydat, 
                control.family = list(hyper = list(theta1=mytheta1, theta2=mytheta2)),E=mydat$ld,
                control.fixed= list(prec.intercept = 1),
                control.inla = list(strategy="gaussian"),
                control.mode = list( result = result_1, restart = TRUE),
                working.directory = paste("~/path/results/",k,sep=""),
                num.threads = 35)


rw_list <- list()
rw_sd_list <- list()
regGR <- GRangesList()
for(p in 1:4)
{
  rw <- list()
  rw_sd <- list()
  for(i in c(1:22,"X"))
  {
    k <- paste("chr",i,sep="")
    result <- inla.collect.results(paste("~/path/10kb/results/",k,sep=""))
    estbin <- matrix(result$summary.random$geneid[,2],ncol=4,byrow=F)
    estbin_low <- matrix(result$summary.random$geneid[,4],ncol=4,byrow=F)
    estbin_high <- matrix(result$summary.random$geneid[,6],ncol=4,byrow=F)
    estbin_sd <- matrix(result$summary.random$geneid[,3],ncol=4,byrow=F)
    rw[[k]] <- estbin[,p]
    rw_sd[[k]] <- estbin_sd[,p]
    
    if(p==1){
      minpos <- 0
      reg_size <- 10000
      maxpos <-nrow(estbin)*reg_size
      nregions <- ceiling((maxpos-minpos) / reg_size)
      regstarts <- seq(minpos,(nregions-1)*reg_size+minpos,reg_size)+1
      regends <-  seq(minpos+reg_size,nregions*reg_size+minpos,reg_size)
      regionGR <- GRanges(seqnames=Rle(rep(k,nregions)),ranges=IRanges(regstarts,regends))
      regionGR$id <- paste(rep(k,nregions),":",regstarts,"..",regends,sep="")
      regionGR <- regionGR[1:nrow(estbin)]  
      regGR[[k]] <- regionGR
    }
  }  
  rw_list[[p]] <- rw
  rw_sd_list[[p]] <- rw_sd
}

unlist(regGR)

regGR <- unlist(regGR)

rw_list_iid <- list()
rw_sd_list_iid <- list()
for(p in 1:4)
{
  rw_iid <- list()
  rw_sd_iid <- list()
  for(i in c(1:22,"X"))
  {
    k <- paste("chr",i,sep="")
    result <- inla.collect.results(paste("~/path/results/",k,sep=""))
    estbin <- matrix(result$summary.random$geneid2[,2],ncol=4,byrow=F)
    estbin_low <- matrix(result$summary.random$geneid2[,4],ncol=4,byrow=F)
    estbin_high <- matrix(result$summary.random$geneid2[,6],ncol=4,byrow=F)
    estbin_sd <- matrix(result$summary.random$geneid2[,3],ncol=4,byrow=F)
    rw_iid[[k]] <- estbin[,p]
    rw_sd_iid[[k]] <- estbin_sd[,p]
  }  
  rw_list_iid[[p]] <- rw_iid
  rw_sd_list_iid[[p]] <- rw_sd_iid
}

colnames(estbin) <- c("B%20lymphoblastoid%20cell%20line%3a%20GM12878%20ENCODE%2c%20biol", 
                      "chronic%20myelogenous%20leukemia%20cell%20line%3aK562%20ENCODE%2c%20biol", 
                      "epitheloid%20carcinoma%20cell%20line%3a%20HelaS3%20ENCODE%2c%20biol", 
                      "hepatocellular%20carcinoma%20cell%20line%3a%20HepG2%20ENCODE%2c%20biol")

rws <- do.call("cbind",lapply(rw_list,unlist))
rws_iid <-  do.call("cbind",lapply(rw_list_iid,unlist))
rws_sd <- do.call("cbind",lapply(rw_sd_list,unlist))
rws_sd_iid <- do.call("cbind",lapply(rw_sd_list_iid,unlist))

colnames(rws) <- c("GLM12878","K562","HelaS3","HepG2")
colnames(rws_iid) <- c("GLM12878","K562","HelaS3","HepG2")
colnames(rws_sd) <- c("GLM12878","K562","HelaS3","HepG2")
colnames(rws_sd_iid) <- c("GLM12878","K562","HelaS3","HepG2")

save(rws,"~/path/results/rws_ENCODE.Rdat")
save(rws_iid,"~/path/results/rws_ENCODE.Rdat")
save(rws_sd,"~/path/results/rws_ENCODE.Rdat")
save(rws_sd_iid,"~/path/results/rws_ENCODE.Rdat")
