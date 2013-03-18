kfoldCV <-
function(times,obs,n,b,K,delta,lambda,annot,annot_node,kfold,active_mu,active_sd,
         inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE){
  kds <- matrix(b,nrow=n,ncol=K)
  # define k-fold groups: stratified
  obs_kfold <- list()
  num <- (dim(obs)[1]*dim(obs)[2])
  le <- ceiling(num/kfold)
  sq_err_all <- edges_all <- vector()
  for(k in 1:times){
	tmps <- sample(seq(1,kfold),kfold)
 	sq_err_tmp <- vector()
	for(j in 1:kfold){
	  tmp <- c(tmps[j:kfold],tmps[1:j-1])
	  obs_order <- order(obs)
	  obs_kfold[[j]] <- matrix(NA,nrow=n,ncol=K)
	  for(i in 1:le){
		if(num>=i){
		  obs_kfold[[j]][obs_order[tmp[1]]] <- obs[obs_order[tmp[1]]]
		}
  # 	  print(tmp)
		obs_order <- obs_order[-tmp]
		tmp <- c(tmp[-1],tmp[1])
	  }
	  obs_kfold[[j]] <- matrix(obs_kfold[[j]],nrow=n,ncol=K)
	}
  #   mean(obs_kfold[[1]],na.rm=T)
	## make crossvalidation
	adja_sum <- adja_num <- matrix(0,ncol=n,nrow=n)
	for(x in 1:kfold){
	  test_ids <- seq(1,kfold)[-x]
	  train_tmp <- vector()
	  for(i in test_ids){
		train_tmp <- rbind(train_tmp,c(obs_kfold[[i]]))
	  }
	  train_data <- rep(NA,dim(train_tmp)[2])
	  for(i in 1:dim(train_tmp)[2]){
		if(!all(is.na(train_tmp[,i]))){
		  train_data[i] <- sum(train_tmp[,i],na.rm=T)
		}
	  }
	  train_data <- matrix(train_data,nrow=n,ncol=K)
	  obs_modified <- train_data
	  ## do ILP
	  res <- doILP(obs_modified,delta,lambda=lambda,b,n,K,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos)
	  adja <- getAdja(res,n)
	  ## calculate statistics on learned edges
	  edges_all <- rbind(edges_all,c(t(adja)))
	  
	  ## calculate mean squared error of predicted and observed
	  predict <- calcPredictionKfoldCV(adja,b,n,K,active_mu,active_sd,inactive_mu,inactive_sd)
	  ids_rem <- which(is.na(obs_modified))
	  sq_err_tmp <- c(sq_err_tmp,((predict[ids_rem]-obs[ids_rem])^2))
	}
	sq_err_all <- rbind(sq_err_all,sq_err_tmp)
  }
  sq_err <- apply(sq_err_all,2,mean,na.rm=T)
  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  edges_all <- edges_all[,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err,na.rm=T)
  return(list(MSE=MSE,edges_all=edges_all))
}
