test.runitKfoldCV_timeSeries <- function() {
	
	n <- 3
	K <- 4
	T_ <- 4

	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)

	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

	obs_mat <- array(NA, c(n,K,T_))

	obs_mat[,,1] <- matrix(c(0.56, 0.56, 0.56, 0.56,
													0.56, 0.56, 0.56, 0.56,
													0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,2] <- matrix(c(0.56, 0.95, 0.95, 0.95,
													0.56, 0.56, 0.56, 0.56,
													0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,3] <- matrix(c(0.56, 0.95, 0.95, 0.95,
													 0.56, 0.56, 0.95, 0.95,
													 0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)
													
	obs_mat[,,4] <- matrix(c(0.56, 0.95, 0.95, 0.95,
													 0.56, 0.56, 0.95, 0.95,
													 0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)

	baseline <- c(0.76,0.76,0)
		
	mu_types <- c("single", "perGene", "perGeneExp", "perGeneTime", "perGeneExpTime")
	delta_types <- c("perGene", "perGene", "perGeneExp", "perGeneTime", "perGeneExpTime")

	mu_list <- list()
	mu_list[[1]] <- list()
	mu_list[[2]] <- list()
	mu_list[[3]] <- list()
	mu_list[[4]] <- list()
	mu_list[[5]] <- list()

	mu_list[[1]]$active_mu <- 0.95
	mu_list[[1]]$active_sd <- 0.01
	mu_list[[1]]$inactive_mu <- 0.56
	mu_list[[1]]$inactive_sd <- 0.01
	mu_list[[1]]$delta <- rep(0.755, n)


	mu_list[[2]]$active_mu <- rep(0.95, n)
	mu_list[[2]]$active_sd <- rep(0.01, n)
	mu_list[[2]]$inactive_mu <- rep(0.56, n)
	mu_list[[2]]$inactive_sd <- rep(0.01, n)
	mu_list[[2]]$delta <- rep(0.755, n)

	mu_list[[3]]$active_mu <- matrix(rep(0.95, n*K), nrow=n, ncol=K)
	mu_list[[3]]$active_sd <- matrix(rep(0.01, n*K), nrow=n, ncol=K)
	mu_list[[3]]$inactive_mu <- matrix(rep(0.56, n*K), nrow=n, ncol=K)
	mu_list[[3]]$inactive_sd <- matrix(rep(0.01, n*K), nrow=n, ncol=K)
	mu_list[[3]]$delta <- matrix(rep(0.755, n*K), nrow=n, ncol=K)

	mu_list[[4]]$active_mu <- matrix(rep(0.95, n*T_), nrow=n, ncol=T_)
	mu_list[[4]]$active_sd <- matrix(rep(0.01, n*T_), nrow=n, ncol=T_)
	mu_list[[4]]$inactive_mu <- matrix(rep(0.56, n*T_), nrow=n, ncol=T_)
	mu_list[[4]]$inactive_sd <- matrix(rep(0.01, n*T_), nrow=n, ncol=T_)
	mu_list[[4]]$delta <- matrix(rep(0.755, n*T_), nrow=n, ncol=T_)

	mu_list[[5]]$active_mu <- array(rep(0.95, n*K*T_), c(n,K,T_))
	mu_list[[5]]$active_sd <- array(rep(0.01, n*K*T_), c(n,K,T_))
	mu_list[[5]]$inactive_mu <- array(rep(0.56, n*K*T_), c(n,K,T_))
	mu_list[[5]]$inactive_sd <- array(rep(0.01, n*K*T_), c(n,K,T_))
	mu_list[[5]]$delta <- array(rep(0.755, n*K*T_), c(n,K,T_))

	kfold <- 10
	lambda <- 1/10
	annot <-  getEdgeAnnot(n)
	annot_node <- seq(1,n)
	
	true_result <- matrix(c(0, 0.7947368, 0.0000000, 
													0, 0.0000000, 0.7947368, 
													0, 0.0000000, 0.0000000), nrow=n, ncol=n, byrow=TRUE)
													
	colnames(true_result) <- rownames(true_result) <-  seq(1,n)

	for (i in 1:length(mu_types)) {
		mu_type <- mu_types[i]
		delta_type <- delta_types[i]
		
		active_mu <- mu_list[[i]]$active_mu
		active_sd <- mu_list[[i]]$active_sd
		inactive_mu <- mu_list[[i]]$inactive_mu
		inactive_sd <- mu_list[[i]]$inactive_sd
		delta <- mu_list[[i]]$delta

		res <- kfoldCV(kfold=kfold, times=1, obs=obs_mat, delta=delta, lambda=lambda, b=b, n=n, K=K, T_=T_, annot=annot,
									 annot_node=annot_node, active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu,
									 inactive_sd=inactive_sd, mu_type=mu_type, delta_type=delta_type, prior=NULL, sourceNode=NULL, 
									 sinkNode=NULL, allint=FALSE, allpos=FALSE, flag_time_series=TRUE)

		adja <- getSampleAdjaMAD(res$edges_all, n, annot_node, method=median, method2=mad, septype="->")
		checkEquals(true_result, adja, tolerance=0.00001)
	}
}

