loocv <-
function(times,obs,n,b,K,delta,lambda,annot,annot_node,
	active_mu,active_sd,inactive_mu,inactive_sd,prior=NULL,sourceNode=NULL,sinkNode=NULL,allint=FALSE,allpos=FALSE){
  kds <- matrix(b,nrow=n,ncol=K)
  # elements to leave out (each element at least once)
  looc <- cbind(rep(seq(1,dim(obs)[1]),dim(obs)[2]),rep(seq(1,dim(obs)[2]),rep(dim(obs)[1],dim(obs)[2])))
  adja_sum <- adja_num <- matrix(0,ncol=n,nrow=n)
  edges_all <- sq_err <- vector()
  # observation of genes n after knockdowns k 
  for(x in 1:dim(looc)[1]){
	sq_err_tmp <- vector()
	obs_modified <- obs
	## randomly select an entry to be missing
	rem_gene <- looc[x,1] 	# which gene is removed
	rem_kd <- looc[x,2]		# in which knockdown
	obs_modified[rem_gene,rem_kd] <- NA
	ele <- obs[rem_gene,rem_kd]
	# mache nur, wenn der datenpunkt nicht eh schon NA ist
	if(!is.na(ele)){
	  ## do ILP
	  res <- doILP(obs_modified,delta,lambda=lambda,b,n,K,annot,prior=prior,sourceNode=sourceNode,sinkNode=sinkNode,all.int=allint,all.pos=allpos)
	  adja <- getAdja(res,n)
	  for(i in 1:times){
		## calculate mean squared error of predicted and observed
		predict <- calcPredictionLOOCV(kds,adja,obs_modified,delta,rem_kd,rem_gene,
		  active_mu,active_sd,inactive_mu,inactive_sd)
		sq_err_tmp <- c(sq_err_tmp,((predict-ele)^2))
	  }
	  res <- NA
	  ## calculate statistics on learned edges
	  edges_all <- rbind(edges_all,c(t(adja)))
	  sq_err <- c(sq_err,mean(sq_err_tmp,na.rm=T))
	}
  }
#   adja_mu <- adja_sum/(times*dim(looc)[1])
#   adja_prob <- adja_num/(times*dim(looc)[1])
  tmp1 <- rep(annot_node,rep(n,n))
  tmp2 <- rep(annot_node,n)
  id_selfloop <- which(tmp1==tmp2)
  tmp <- paste(tmp1,tmp2,sep="->")
  edges_all <- edges_all[,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err,na.rm=T)
  return(list(MSE=MSE,edges_all=edges_all))
}
