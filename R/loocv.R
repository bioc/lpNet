#
# do leave one out cross validation
#
loocv <- function(kfold=NULL, times, obs, delta, lambda, b, n, K, T_=NULL, annot, annot_node, active_mu, active_sd, 
									inactive_mu, inactive_sd, mu_type, delta_type, prior=NULL, sourceNode=NULL, sinkNode=NULL, 
									allint=FALSE, allpos=FALSE, flag_time_series=FALSE) {
										
	if (flag_time_series == FALSE) {
		res <- .loocv_steadyState(kfold, times, obs, delta, lambda, b, n, K, T_, annot, annot_node, 
															active_mu, active_sd, inactive_mu, inactive_sd, mu_type, delta_type, 
															prior, sourceNode, sinkNode, allint, allpos)
	}
	else {
		res <- .loocv_timeSeries(kfold, times, obs, delta, lambda, b, n, K, T_, annot, annot_node,
														 active_mu, active_sd, inactive_mu, inactive_sd, mu_type, delta_type,
														 prior, sourceNode, sinkNode, allint, allpos)
	}
	
	return(res)
}


.loocv_steadyState <- function(kfold=NULL, times, obs, delta, lambda, b, n, K, T_=NULL, annot, annot_node,
															 active_mu, active_sd, inactive_mu, inactive_sd, mu_type, delta_type, 
															 prior=NULL, sourceNode=NULL, sinkNode=NULL, allint=FALSE, allpos=FALSE) {

  # elements to leave out (each element at least once)
  looc <- cbind(rep(seq(1, dim(obs)[1]), dim(obs)[2]),
								rep(seq(1, dim(obs)[2]), rep(dim(obs)[1], dim(obs)[2])))

  edges_all <- sq_err <- baseline_all <- vector()
  
  # observation of genes n after knockdowns k 
  for (x in 1:dim(looc)[1]) {
		sq_err_tmp <- vector()
		obs_modified <- obs

		rem_gene <- looc[x,1] 	# which gene is removed
		rem_kd <- looc[x,2]		# in which knockdown
		obs_modified[rem_gene,rem_kd] <- NA
		ele <- obs[rem_gene,rem_kd]

		if (!is.na(ele)) {
			## do ILP
			res <- .doILP_steadyState(obs=obs_modified, delta=delta, lambda=lambda, b=b, n=n, K=K, T_=T_, annot=annot, 
																delta_type=delta_type,prior=prior, sourceNode=sourceNode, sinkNode=sinkNode, 
																all.int=allint, all.pos=allpos)
						
			adja <- getAdja(res=res, n=n)
			baseline <- getBaseline(res=res, n=n)
			edges_all <- rbind(edges_all, c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)

			for (i in 1:times) {  # calculate prediction a given number of times
				predict <- .calcPredictionLOOCV_steadyState(obs=obs_modified, delta=delta, b=b, n=n, K=K, adja=adja, baseline=baseline, 
																										rem_gene=rem_gene, rem_k=rem_kd, active_mu=active_mu, 
																										active_sd=active_sd, inactive_mu=inactive_mu, inactive_sd=inactive_sd, 
																										mu_type=mu_type)
				
				sq_err_tmp <- c(sq_err_tmp, ((predict-ele)^2))  # calculate mean squared error of predicted and observed
			}
			sq_err <- c(sq_err,mean(sq_err_tmp,na.rm=T))  # calculate statistics on learned edges
		}
  }


  tmp1 <- rep(annot_node, rep(n, n))
  tmp2 <- rep(annot_node, n)
  id_selfloop <- which(tmp1 == tmp2)
  tmp <- paste(tmp1, tmp2, sep="->")
  edges_all <- edges_all[, -id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err, na.rm=T)
  
  return(list(MSE=MSE, edges_all=edges_all, baseline_all=baseline_all))
}


.loocv_timeSeries <-function(kfold=NULL, times, obs, delta, lambda, b, n, K, T_, annot, annot_node, 
														 active_mu, active_sd, inactive_mu, inactive_sd, mu_type, delta_type,
														 prior=NULL, sourceNode=NULL, sinkNode=NULL, allint=FALSE, allpos=FALSE) {
	
  # elements to leave out (each element at least once)
  looc <- cbind(rep(seq(1, dim(obs)[1]), dim(obs)[2]),
								rep(seq(1, dim(obs)[2]), rep(dim(obs)[1], dim(obs)[2])),
								rep(seq(1, dim(obs)[3]), rep(dim(obs)[1] * dim(obs)[2], dim(obs)[3])))
								
  edges_all <- sq_err <- baseline_all <- vector()
	
	# dont remove observations from first time point
  looc = looc[-which(looc[ ,3] == 1, arr.ind=T),]
	
	# observation of genes n after knockdowns k 
	for (x in 1:dim(looc)[1]) {
		sq_err_tmp <- vector()
		obs_modified <- obs
		rem_gene <- looc[x,1]  # which gene is removed
		rem_kd <- looc[x,2]  # in which knockdown
		rem_t <- looc[x,3]  # in which time point
		obs_modified[rem_gene,rem_kd,rem_t] <- NA
		ele <- obs[rem_gene,rem_kd,rem_t]
		
		if (!is.na(ele)) {
			# do ILP
			res <- .doILP_timeSeries(obs=obs_modified, delta=delta, lambda=lambda, b=b, n=n, K=K, T_=T_, annot=annot, delta_type=delta_type,
															 prior=prior, sourceNode=sourceNode, sinkNode=sinkNode, all.int=allint, all.pos=allpos)
			
			adja <- getAdja(res=res, n=n)
			baseline <- getBaseline(res=res, n=n)
			edges_all <- rbind(edges_all, c(t(adja)))
			baseline_all <- rbind(baseline_all, baseline)
			
			for (i in 1:times) {  # calculate prediction a given number of times
				predict <- .calcPredictionLOOCV_timeSeries(obs=obs_modified, delta=delta, b=b, n=n, K=K, adja=adja, baseline=baseline, 
																									 rem_gene=rem_gene, rem_k=rem_kd, rem_t=rem_t, active_mu=active_mu, 
																									 active_sd=active_sd, inactive_mu=inactive_mu, inactive_sd=inactive_sd, mu_type=mu_type)
				
				sq_err_tmp <- c(sq_err_tmp, ((predict-ele)^2))  # calculate mean squared error of predicted and observed
			}		
			sq_err <- c(sq_err, mean(sq_err_tmp,na.rm=T))  # calculate statistics on learned edges
		}
	}

  tmp1 <- rep(annot_node, rep(n, n))
  tmp2 <- rep(annot_node, n)
  id_selfloop <- which(tmp1 == tmp2)
  tmp <- paste(tmp1, tmp2, sep="->")
  edges_all <- edges_all[ ,-id_selfloop]
  colnames(edges_all) <- tmp[-id_selfloop]
  MSE <- mean(sq_err, na.rm=T)

  return(list(MSE=MSE, edges_all=edges_all, baseline_all=baseline_all))
}
