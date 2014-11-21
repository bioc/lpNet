.setUp <- function() {

	n <<- 3
	K <<- 4

	T_nw <<- matrix(c(0,1,0,
										0,0,1,
										0,0,0), nrow=n, ncol=n, byrow=TRUE)

	b <<- c(0,1,1,
					1,0,1,
					1,1,0,
					1,1,1)

	obs_mat <<- matrix(c(0.56, 0.95, 0.95, 0.95,
											 0.56, 0.56, 0.95, 0.95,
											 0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)

	baseline <<- c(0.76, 0.76, 0)
	
	mu_types <<- c("single", "perGene", "perGeneExp")

	mu_list <<- list()
	mu_list[[1]] <<- list()
	mu_list[[2]] <<- list()
	mu_list[[3]] <<- list()

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
}


test.runitCalcPredictionKfoldCV <- function() {

	obs_modified <- obs_mat
	obs_modified[2,4] <- NA

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
		predict <- calcPredictionKfoldCV(obs, delta, b, n, K, adja=T_nw, baseline, rem_entries, rem_entries_vec,
																		 active_mu, active_sd, inactive_mu, inactive_sd, mu_type=mu_type) 
		
		checkEquals(obs_mat, predict, tolerance=0.05)
	}
}
