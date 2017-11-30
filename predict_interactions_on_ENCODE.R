#Script for running the random forest interactions modelling
#requires test dataset generated from previous scripts

library(GenomicRanges)
library(foreach)
library(doParallel)
library(randomForest)
library(pROC)
library(caret)
library(unbalanced)

registerDoParallel(cores=30)
load(file="/binf-isilon/alab/projects/ARCHS/ARCHS_RW/interactions_modelling/ENCODE_interactions_testdats.Rdat")

myDistRange <- 6:200
glob_use <- unique(unlist(lapply(enc_dat,function(x) which( (x$bait_in_group>0) & (x$dist %in% myDistRange)  ))))

glob_use <- sort(glob_use)
length(glob_use)

mydatarun <- enc_dat[[1]][glob_use,]
mydatarun$sig <- mydatarun$count>=3

enc_all <- lapply(enc_dat,function(x) x[glob_use,])

#mydatarun <- mydatarun[en_supported,]

data <- mydatarun
data$sig <- mydatarun$sig


predictions_lists <- list()
for(myrun in c(1:10))
{
  data <- mydatarun
  data$sig <- mydatarun$sig
  
  pred_list <- list()
  test_sets <- list()
  training_sets <- list()
  
  trainingset <- data
  mysig <- data$sig
  trainingset <- trainingset[,-which(colnames(trainingset) %in% c("bait_in_group","bait_en_support","bait_nens","bait_eRNAsum","bait_CTCF","targ_CTCF","chr","votes_random","votes_top","votes_bottom","probability","sig2","hicbound","raobound","eigcross","count","count","is.min","bait","targ","en_bet","bait_id","targ_id","sig","votes"))]
  trainingset$sig <- mysig
  
  #   ###OVER ALL DISTANCES
  myub <- list()
  
  registerDoParallel(cores=20)
  
  myub <- foreach(d=(min(myDistRange):max(myDistRange))) %dopar% {
    print(d)
    myrfdat <- trainingset[trainingset$dist==d,]
    mysig <- myrfdat$sig
    myrfdat <- myrfdat[,-which(colnames(trainingset) %in% "sig")]
    if( length(mysig[mysig])>0 )
    {
      ubBalance(X=myrfdat,Y=factor(as.numeric(mysig)),positive="1",percOver=200,percUnder=150,k=min(max(sum(mysig),2),5))
    }
  }
  
  myrfdat <- do.call("rbind",lapply(myub,function(x) x$X))
  mysig <- as.numeric(unlist(lapply(myub,function(x) x$Y)))>1
  
  trainingset <- myrfdat
  trainingset$sig <- mysig
  
  registerDoParallel(cores=30)
  
  train.rf <- trainingset[,-which(colnames(trainingset)=="sig")]
  train.sig.rf <- factor(trainingset$sig)
  names(train.sig.rf) <- row.names(train.rf)
  cf1 <- foreach(ntree=rep(15, 30), .combine=combine, .packages='randomForest') %dopar% randomForest(train.rf,train.sig.rf ,norm.votes=TRUE, importance=TRUE, ntree=ntree)  
  
  
  registerDoParallel(cores=5)
  
  pred_list <- foreach(i=(1:4)) %dopar% {
    testset <- enc_dat[[i]][glob_use,]
    print(paste(i,":perform predictions on test data..."))
    predict(cf1, newdata=testset,type="prob",norm.votes=T)[,2]
  }
  
  predictions_lists[[myrun]] <- pred_list
}

predictions_averages <- list()
for(c in 1:4)
{
  print(c)
  predictions_averages[[c]] <- rep(0,length(predictions_lists[[1]][[c]]))
  for(k in 1:10)
  {
    predictions_averages[[c]] <- predictions_averages[[c]] + predictions_lists[[k]][[c]]
  }
}

for(c in 1:4)
{
  predictions_averages[[c]] <- predictions_averages[[c]]/10
}


load(file=paste(serve,"/binf-isilon/alab/projects/ARCHS/ARCHS_RW/interactions_modelling/ENCODE_datasets_with_predicted_interaction_probabilities.Rdat",sep=""))
gm_preds <- pred_sets[[1]]$probability

pred_sets <- list()
for(i in 1:4)
{
  print(i)
  testset <- enc_dat[[i]][glob_use,]
  pred_sets[[i]] <- testset
  pred_sets[[i]]$probability <- predictions_averages[[i]]
  if(i==1)
  {
    pred_sets[[i]]$probability <- gm_preds
  }
}


  save(pred_sets,file="/binf-isilon/alab/projects/ARCHS/ARCHS_RW/interactions_modelling/ENCODE_datasets_with_predicted_interaction_probabilities_10_runs.Rdat")
  
    
  load(file=paste(serve,"/binf-isilon/alab/projects/ARCHS/ARCHS_RW/interactions_modelling/ENCODE_datasets_with_predicted_interaction_probabilities_10_runs.Rdat",sep=""))
  
  en_supported <- which(pred_sets[[1]]$targ_nens>0 | pred_sets[[2]]$targ_nens>0 | pred_sets[[3]]$targ_nens>0 | pred_sets[[4]]$targ_nens>0)
  #lapply(pred_sets,nrow)
  
  
  EP_sets <- list()
  for(i in 1:4)
  {
    EP_sets[[i]] <- pred_sets[[i]][en_supported,]
  }
  lapply(EP_sets,nrow)
  