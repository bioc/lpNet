.setUp <- function() {

	n <<- 3
	K <<- 4
	T_ <<- 3

	T_nw <<- matrix(c(0,0,1,
									  0,0,-1,
									  0,0,0), nrow=n, ncol=n, byrow=TRUE)

	b <<- c(0,1,1,
					1,0,1,
					1,1,0,
					1,1,1)

	obs_mat <<- array(NA, c(n,K,T_))

	obs_mat[,,1] <<- matrix(c(0.56, 0.56, 0.56, 0.56,
														0.56, 0.56, 0.56, 0.56,
														0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,2] <<- matrix(c(0.56, 0.95, 0.95, 0.95,
														0.95, 0.56, 0.95, 0.95,
														0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,3] <<- matrix(c(0.56, 0.95, 0.95, 0.95,
														0.95, 0.56, 0.95, 0.95,
														0.56, 0.95, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)

	baseline <<- c(0.76, 0.76, 0)

	mu_types <<- c("single", "perGene", "perGeneExp", "perGeneTime", "perGeneExpTime")

	mu_list <<- list()
	mu_list[[1]] <<- list()
	mu_list[[2]] <<- list()
	mu_list[[3]] <<- list()
	mu_list[[4]] <<- list()
	mu_list[[5]] <<- list()

	mu_list[[1]]$active_mu <<- 0.95
	mu_list[[1]]$active_sd <<- 0.01
	mu_list[[1]]$inactive_mu <<- 0.56
	mu_list[[1]]$inactive_sd <<- 0.01
	mu_list[[1]]$delta <<- rep(0.755, n)

	mu_list[[2]]$active_mu <<- rep(0.95, n)
	mu_list[[2]]$active_sd <<- rep(0.01, n)
	mu_list[[2]]$inactive_mu <<- rep(0.56, n)
	mu_list[[2]]$inactive_sd <<- rep(0.01, n)
	mu_list[[2]]$delta <<- rep(0.755, n)

	mu_list[[3]]$active_mu <<- matrix(rep(0.95, n*K), nrow=n, ncol=K)
	mu_list[[3]]$active_sd <<- matrix(rep(0.01, n*K), nrow=n, ncol=K)
	mu_list[[3]]$inactive_mu <<- matrix(rep(0.56, n*K), nrow=n, ncol=K)
	mu_list[[3]]$inactive_sd <<- matrix(rep(0.01, n*K), nrow=n, ncol=K)
	mu_list[[3]]$delta <<- matrix(rep(0.755, n*K), nrow=n, ncol=K)

	mu_list[[4]]$active_mu <<- matrix(rep(0.95, n*T_), nrow=n, ncol=T_)
	mu_list[[4]]$active_sd <<- matrix(rep(0.01, n*T_), nrow=n, ncol=T_)
	mu_list[[4]]$inactive_mu <<- matrix(rep(0.56, n*T_), nrow=n, ncol=T_)
	mu_list[[4]]$inactive_sd <<- matrix(rep(0.01, n*T_), nrow=n, ncol=T_)
	mu_list[[4]]$delta <<- matrix(rep(0.755, n*T_), nrow=n, ncol=T_)

	mu_list[[5]]$active_mu <<- array(rep(0.95, n*K*T_), c(n,K,T_))
	mu_list[[5]]$active_sd <<- array(rep(0.01, n*K*T_), c(n,K,T_))
	mu_list[[5]]$inactive_mu <<- array(rep(0.56, n*K*T_), c(n,K,T_))
	mu_list[[5]]$inactive_sd <<- array(rep(0.01, n*K*T_), c(n,K,T_))
	mu_list[[5]]$delta <<- array(rep(0.755, n*K*T_), c(n,K,T_))
}


test.runitCalcPredictionLOOCV01 <- function() {
	
	T_nw <- matrix(c(0,0,1,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
	
	obs_mat[,,1] <- matrix(c(0.56, 0.56, 0.56, 0.56,
													 0.56, 0.56, 0.56, 0.56,
													 0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,2] <- matrix(c(0.56, 0.95, 0.95, 0.95,
													 0.95, 0.56, 0.95, 0.95,
													 0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,3] <- matrix(c(0.56, 0.95, 0.95, 0.95,
													 0.95, 0.56, 0.95, 0.95,
													 0.95, 0.95, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)
	baseline <- c(0, 0, 0)
	
	obs_modified <- obs_mat
	rem_gene <- 2
	rem_k <- 4
	rem_t <- 2
	obs_modified[2,4,2] <- NA

	rem_entries <- which(is.na(obs_modified), arr.ind=TRUE)
	rem_entries_vec <- which(is.na(obs_modified))
	
	for (i in 1:length(mu_types)) {
		mu_type <- mu_types[i]
		active_mu <- mu_list[[i]]$active_mu
		active_sd <- mu_list[[i]]$active_sd
		inactive_mu <- mu_list[[i]]$inactive_mu
		inactive_sd <- mu_list[[i]]$inactive_sd
		delta <- mu_list[[i]]$delta
		
		## calculate mean squared error of predicted and observed
		predict <- calcPredictionLOOCV(obs=obs_modified, delta=delta, b=b, n=n, K=K, adja=T_nw, 
																	 baseline=baseline, rem_gene=rem_gene, rem_k=rem_k, rem_t=rem_t,
																	 active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu, 
																	 inactive_sd=inactive_sd, mu_type=mu_type, flag_time_series=TRUE)
		
		checkEquals(predict, 0.56, tolerance=0.05)
	}
}


test.runitCalcPredictionLOOCV02 <- function() {
	
	T_nw <- matrix(c(0,0,1,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
									 
	obs_mat[,,1] <- matrix(c(0.56, 0.56, 0.56, 0.56,
													 0.56, 0.56, 0.56, 0.56,
													 0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,2] <- matrix(c(0.56, 0.95, 0.95, 0.95,
													 0.95, 0.56, 0.95, 0.95,
													 0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,3] <- matrix(c(0.56, 0.95, 0.95, 0.95,
													 0.95, 0.56, 0.95, 0.95,
													 0.95, 0.95, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)
	
	obs_modified <- obs_mat
	rem_gene <- 2
	rem_k <- 4
	rem_t <- 2
	obs_modified[2,4,2] <- NA

	rem_entries <- which(is.na(obs_modified), arr.ind=TRUE)
	rem_entries_vec <- which(is.na(obs_modified))
	
	
	for (i in 1:length(mu_types)) {
		mu_type <- mu_types[i]
		active_mu <- mu_list[[i]]$active_mu
		active_sd <- mu_list[[i]]$active_sd
		inactive_mu <- mu_list[[i]]$inactive_mu
		inactive_sd <- mu_list[[i]]$inactive_sd
		delta <- mu_list[[i]]$delta
		predict <- calcPredictionLOOCV(obs=obs_modified, delta=delta, b=b, n=n, K=K, adja=T_nw, 
																	 baseline=baseline, rem_gene=rem_gene, rem_k=rem_k, rem_t=rem_t,
																	 active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu,
																	 inactive_sd=inactive_sd, mu_type=mu_type, flag_time_series=TRUE)
		
		checkEquals(predict, 0.95, tolerance=0.05)
	}
}


test.runitCalcPredictionLOOCV03 <- function() {
	
	T_nw <- matrix(c(0,0,1,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)

	obs_mat[,,1] <- matrix(c(0.56, 0.56, 0.56, 0.56,
													 0.56, 0.56, 0.56, 0.56,
													 0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,2] <- matrix(c(0.56, 0.95, 0.95, 0.95,
													 0.95, 0.56, 0.95, 0.95,
													 0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,3] <- matrix(c(0.56, 0.95, 0.95, 0.95,
													 0.95, 0.56, 0.95, 0.95,
													 0.95, 0.95, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)
	
	obs_modified <- obs_mat
	rem_gene <- 3
	rem_k <- 4
	rem_t <- 3
	obs_modified[3,4,3] <- NA

	rem_entries <- which(is.na(obs_modified), arr.ind=TRUE)
	rem_entries_vec <- which(is.na(obs_modified))
	
	for (i in 1:length(mu_types)) {
		mu_type <- mu_types[i]
		active_mu <- mu_list[[i]]$active_mu
		active_sd <- mu_list[[i]]$active_sd
		inactive_mu <- mu_list[[i]]$inactive_mu
		inactive_sd <- mu_list[[i]]$inactive_sd
		delta <- mu_list[[i]]$delta
		
		predict <- calcPredictionLOOCV(obs=obs_modified, delta=delta, b=b, n=n, K=K, adja=T_nw, 
																	 baseline=baseline, rem_gene=rem_gene, rem_k=rem_k, rem_t=rem_t,
																	 active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu, 
																	 inactive_sd=inactive_sd, mu_type=mu_type, flag_time_series=TRUE)

	checkEquals(predict, 0.95, tolerance=0.05)
	}
}


test.runitCalcPredictionLOOCV04 <- function() {
	
	T_nw <- matrix(c(0,0,1,
									 0,0,-1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
		
	obs_modified <- obs_mat
	rem_gene <- 3
	rem_k <- 4
	rem_t <- 3
	obs_modified[2,4,2] <- NA
	obs_modified[3,4,3] <- NA

	rem_entries <- which(is.na(obs_modified), arr.ind=TRUE)
	rem_entries_vec <- which(is.na(obs_modified))
	
	
	for (i in 1:length(mu_types)) {
		mu_type <- mu_types[i]
		active_mu <- mu_list[[i]]$active_mu
		active_sd <- mu_list[[i]]$active_sd
		inactive_mu <- mu_list[[i]]$inactive_mu
		inactive_sd <- mu_list[[i]]$inactive_sd
		delta <- mu_list[[i]]$delta
		
		predict <- calcPredictionLOOCV(obs=obs_modified, delta=delta,  b=b, n=n, K=K, adja=T_nw, baseline=baseline, 
																	 rem_gene=rem_gene, rem_k=rem_k, rem_t=rem_t,
																	 active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu, 
																	 inactive_sd=inactive_sd, mu_type=mu_type, flag_time_series=TRUE)

		checkTrue(is.na(predict))
	}
}


test.runitCalcPredictionLOOCV05 <- function() {

	obs_modified <- obs_mat
	rem_gene <- 3
	rem_k <- 2
	rem_t <- 3
	obs_modified[2,2,2] <- NA
	obs_modified[3,2,3] <- NA

	rem_entries <- which(is.na(obs_modified), arr.ind=TRUE)
	rem_entries_vec <- which(is.na(obs_modified))
	
	for (i in 1:length(mu_types)) {
		mu_type <- mu_types[i]
		active_mu <- mu_list[[i]]$active_mu
		active_sd <- mu_list[[i]]$active_sd
		inactive_mu <- mu_list[[i]]$inactive_mu
		inactive_sd <- mu_list[[i]]$inactive_sd
		delta <- mu_list[[i]]$delta
		
		predict <- calcPredictionLOOCV(obs=obs_modified, delta=delta,  b=b, n=n, K=K, adja=T_nw, baseline=baseline, 
																	 rem_gene=rem_gene, rem_k=rem_k, rem_t=rem_t,
																	 active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu,
																	 inactive_sd=inactive_sd, mu_type=mu_type, flag_time_series=TRUE)

	checkEquals(predict, 0.95, tolerance=0.05)
	}
}


test.runitCalcPredictionLOOCV06 <- function() {

	obs_modified <- obs_mat
	rem_gene <- 3
	rem_k <- 2
	rem_t <- 2
	obs_modified[2,2,1] <- NA
	obs_modified[3,2,2] <- NA

	rem_entries <- which(is.na(obs_modified), arr.ind=TRUE)
	rem_entries_vec <- which(is.na(obs_modified))

	for (i in 1:length(mu_types)) {
		mu_type <- mu_types[i]
		active_mu <- mu_list[[i]]$active_mu
		active_sd <- mu_list[[i]]$active_sd
		inactive_mu <- mu_list[[i]]$inactive_mu
		inactive_sd <- mu_list[[i]]$inactive_sd
		delta <- mu_list[[i]]$delta
		
		predict <- calcPredictionLOOCV(obs=obs_modified, delta=delta, b=b, n=n, K=K, adja=T_nw, 
																	 baseline=baseline, rem_gene=rem_gene, rem_k=rem_k, rem_t=rem_t,
																	 active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu, 
																	 inactive_sd=inactive_sd, mu_type=mu_type, flag_time_series=TRUE)

		checkTrue(is.na(predict))
	}
}


test.runitCalcPredictionLOOCV07 <- function() {

	baseline <- c(0.76, 0.76, 0.76)

	obs_modified <- obs_mat
	rem_gene <- 3
	rem_k <- 2
	rem_t <- 2
	obs_modified[2,2,1] <- NA
	obs_modified[3,2,2] <- NA

	rem_entries <- which(is.na(obs_modified), arr.ind=TRUE)
	rem_entries_vec <- which(is.na(obs_modified))

	for (i in 1:length(mu_types)) {
		mu_type <- mu_types[i]
		active_mu <- mu_list[[i]]$active_mu
		active_sd <- mu_list[[i]]$active_sd
		inactive_mu <- mu_list[[i]]$inactive_mu
		inactive_sd <- mu_list[[i]]$inactive_sd
		delta <- mu_list[[i]]$delta
		
		predict <- calcPredictionLOOCV(obs=obs_modified, delta=delta, b=b, n=n, K=K, adja=T_nw, 
																	 baseline=baseline, rem_gene=rem_gene, rem_k=rem_k, rem_t=rem_t,
																	 active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu, 
																	 inactive_sd=inactive_sd, mu_type=mu_type, flag_time_series=TRUE)
																								
		checkEquals(predict, 0.95, tolerance=0.05)
	}
}
