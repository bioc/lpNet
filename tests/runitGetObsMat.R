test.getObsMatMuTypeSingle <- function() {

	n <- 3
	K <- 4
	
	true_result <- matrix(c(0.56, 0.95, 0.95, 0.95,
													0.56, 0.56, 0.95, 0.95,
													0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=T)

	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

	act_mat <- calcActivation(T_nw, b, n, K)
	
	active_mu <- 0.95
	active_sd <- 0.01
	inactive_mu <- 0.56
	inactive_sd <- 0.01
	
	obs_mat <- getObsMat(act_mat, net_states=NULL, active_mu, active_sd, inactive_mu, inactive_sd, mu_type="single")
	checkEquals(true_result, obs_mat, tolerance=(active_sd + inactive_sd))
}


test.getObsMatMuTypePerGene <- function() {

	n <- 3
	K <- 4
	
	true_result <- matrix(c(0.56, 0.95, 0.95, 0.95,
													0.4, 0.4, 1.1, 1.1,
													0.2, 0.2, 0.2, 1.3), nrow=n, ncol=K, byrow=T)
	
	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

	act_mat <- calcActivation(T_nw, b, n, K)
	
	
	active_mu <- c(0.95, 1.1, 1.3)
	active_sd <- rep(0.01, n)
	inactive_mu <- c(0.56, 0.4, 0.2)
	inactive_sd <- rep(0.01, n)
	
	obs_mat <- getObsMat(act_mat, net_states=NULL, active_mu, active_sd, inactive_mu, inactive_sd, mu_type="perGene")
	checkEquals(true_result, obs_mat, tolerance=(max(active_sd) + max(inactive_sd)))
}


test.getObsMatMuTypePerGeneExp <- function() {

	n <- 3
	K <- 4
	
	true_result <- matrix(c(1.1, 10.3, 10.5, 10.7,
													2.1, 2.3, 20.5, 20.7,
													3.1, 3.3, 3.5, 30.7), nrow=n, ncol=K, byrow=T)
	
	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

	act_mat <- calcActivation(T_nw, b, n, K)
	
	active_mu <- matrix(c(10.1, 20.1, 30.1,
												10.3, 20.3, 30.3,
												10.5, 20.5, 30.5,
												10.7, 20.7, 30.7), nrow=n, ncol=K)
												
	active_sd <- matrix(rep(0.01, n*K), nrow=n, ncol=K)
	
	inactive_mu <- matrix(c(1.1, 2.1, 3.1,
													1.3, 2.3, 3.3,
													1.5, 2.5, 3.5,
													1.7, 2.7, 3.7), nrow=n, ncol=K)
													
	inactive_sd <- matrix(rep(0.01, n*K), nrow=n, ncol=K)
	
	obs_mat <- getObsMat(act_mat, net_states=NULL, active_mu, active_sd, inactive_mu, inactive_sd, mu_type="perGeneExp")
	checkEquals(true_result, obs_mat, tolerance=(max(active_sd) + max(inactive_sd)))
}


test.getObsMatMuTypeSingle_nodeStates <- function() {

	n <- 3
	K <- 4
  T_ <- 4
	
	true_result <- array(NA, c(n, K, T_))
    
	true_result[,,1] <- matrix(c(0.56, 0.56, 0.56, 0.56,
															 0.56, 0.56, 0.56, 0.56,
															 0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=T)
															
	true_result[,,2] <- matrix(c(0.56, 0.95, 0.95, 0.95,
															 0.56, 0.56, 0.56, 0.56,
															 0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=T)

	true_result[,,3] <- matrix(c(0.56, 0.95, 0.95, 0.95,
															 0.56, 0.56, 0.95, 0.95,
															 0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=T)

	true_result[,,4] <- matrix(c(0.56, 0.95, 0.95, 0.95,
															 0.56, 0.56, 0.95, 0.95,
															 0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=T)
	
	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

    net_states <- array(NA, c(n,K,T_))
    
    net_states[,,1] <- matrix(c(0,0,0,0,
                                0,0,0,0,
                                0,0,0,0), nrow=n, ncol=K, byrow=T)
    
	net_states[,,2] <- matrix(c(0,1,1,1,
                                0,0,0,0,
                                0,0,0,0), nrow=n, ncol=K, byrow=T)
	
    net_states[,,3] <- matrix(c(0,1,1,1,
                                0,0,1,1,
                                0,0,0,0), nrow=n, ncol=K, byrow=T)
    
    net_states[,,4] <- matrix(c(0,1,1,1,
                                0,0,1,1,
                                0,0,0,1), nrow=n, ncol=K, byrow=T)
                                
	active_mu <- 0.95
	active_sd <- 0.01
	inactive_mu <- 0.56
	inactive_sd <- 0.01
	
	obs_mat <- getObsMat(act_mat=NULL, net_states, active_mu, active_sd, inactive_mu, inactive_sd, mu_type="single")
	checkEquals(true_result, obs_mat, tolerance=(active_sd + inactive_sd))
}


test.getObsMatMuTypePerGene_nodeStates <- function() {

	n <- 3
	K <- 4
	T_ <- 4
	
	true_result <- array(NA, c(n,K,T_))
	
	true_result[,,1] <- matrix(c(0.56, 0.56, 0.56, 0.56,
															0.4, 0.4, 0.4, 0.4,
															0.2, 0.2, 0.2, 0.2), nrow=n, ncol=K, byrow=T)
	
	true_result[,,2] <- matrix(c(0.56, 0.95, 0.95, 0.95,
															0.4, 0.4, 0.4, 0.4,
															0.2, 0.2, 0.2, 0.2), nrow=n, ncol=K, byrow=T)
	
	true_result[,,3] <- matrix(c(0.56, 0.95, 0.95, 0.95,
															0.4, 0.4, 1.1, 1.1,
															0.2, 0.2, 0.2, 0.2), nrow=n, ncol=K, byrow=T)
															
	true_result[,,4] <- matrix(c(0.56, 0.95, 0.95, 0.95,
															0.4, 0.4, 1.1, 1.1,
															0.2, 0.2, 0.2, 1.3), nrow=n, ncol=K, byrow=T)
	
    
	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)
                 
	net_states <- array(NA, c(n,K,T_))
	
	net_states[,,1] <- matrix(c(0,0,0,0,
															0,0,0,0,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
    
	net_states[,,2] <- matrix(c(0,1,1,1,
															0,0,0,0,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
	
	net_states[,,3] <- matrix(c(0,1,1,1,
															0,0,1,1,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
	
	net_states[,,4] <- matrix(c(0,1,1,1,
															0,0,1,1,
															0,0,0,1), nrow=n, ncol=K, byrow=T)
	
	active_mu <- c(0.95, 1.1, 1.3)
	active_sd <- rep(0.01, n)
	inactive_mu <- c(0.56, 0.4, 0.2)
	inactive_sd <- rep(0.01, n)
	
	obs_mat <- getObsMat(act_mat=NULL, net_states, active_mu, active_sd, inactive_mu, inactive_sd, mu_type="perGene")
	checkEquals(true_result, obs_mat, tolerance=(max(active_sd) + max(inactive_sd)))
}


test.getObsMatMuTypePerGeneExp_nodeStates <- function() {

	n <- 3
	K <- 4
	T_ <- 4
	
	true_result <- array(NA, c(n,K,T_))
    
	true_result[,,1] <- matrix(c(1.1, 1.3, 1.5, 1.7,
															2.1, 2.3, 2.5, 2.7,
															3.1, 3.3, 3.5, 3.7), nrow=n, ncol=K, byrow=T)
															
	true_result[,,2] <- matrix(c(1.1, 10.3, 10.5, 10.7,
															2.1, 2.3, 2.5, 2.7,
															3.1, 3.3, 3.5, 3.7), nrow=n, ncol=K, byrow=T)
															
	true_result[,,3] <- matrix(c(1.1, 10.3, 10.5, 10.7,
															2.1, 2.3, 20.5, 20.7,
															3.1, 3.3, 3.5, 3.7), nrow=n, ncol=K, byrow=T)
															
	true_result[,,4] <- matrix(c(1.1, 10.3, 10.5, 10.7,
															2.1, 2.3, 20.5, 20.7,
															3.1, 3.3, 3.5, 30.7), nrow=n, ncol=K, byrow=T)

	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
									 
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

	net_states <- array(NA, c(n,K,T_))
    
	net_states[,,1] <- matrix(c(0,0,0,0,
															0,0,0,0,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
    
	net_states[,,2] <- matrix(c(0,1,1,1,
															0,0,0,0,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
	
	net_states[,,3] <- matrix(c(0,1,1,1,
															0,0,1,1,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
	
	net_states[,,4] <- matrix(c(0,1,1,1,
															0,0,1,1,
															0,0,0,1), nrow=n, ncol=K, byrow=T)
	
	
	active_mu <- matrix(c(10.1, 10.3, 10.5, 10.7,
                         20.1, 20.3, 20.5, 20.7,
                         30.1, 30.3, 30.5, 30.7), nrow=n, ncol=K, byrow=T)
                         
	active_sd <- matrix(rep(0.01, n*K), nrow=n, ncol=K)
	
	inactive_mu <- matrix(c(1.1, 1.3, 1.5, 1.7,
													 2.1, 2.3, 2.5, 2.7,
													 3.1, 3.3, 3.5, 3.7), nrow=n, ncol=K, byrow=T)

	inactive_sd <- matrix(rep(0.01, n*K), nrow=n, ncol=K)
	
	obs_mat <- getObsMat(act_mat=NULL, net_states,  active_mu, active_sd, inactive_mu, inactive_sd, mu_type="perGeneExp")
	checkEquals(true_result, obs_mat, tolerance=(max(active_sd) + max(inactive_sd)))
}


test.getObsMatMuTypePerGeneTime_nodeStates <- function() {

	n <- 3
	K <- 4
	T_ <- 4
	
	true_result <- array(NA, c(n,K,T_))
    
	true_result[,,1] <- matrix(c(1.1, 1.1, 1.1, 1.1,
															2.1, 2.1, 2.1, 2.1,
															3.1, 3.1, 3.1, 3.1), nrow=n, ncol=K, byrow=T)
															
	true_result[,,2] <- matrix(c(1.3, 10.3, 10.3, 10.3,
															2.1, 2.3, 2.3, 2.3,
															3.3, 3.3, 3.3, 3.3), nrow=n, ncol=K, byrow=T)
															
	true_result[,,3] <- matrix(c(1.5, 10.5, 10.5, 10.5,
															2.5, 2.5, 20.5, 20.5,
															3.5, 3.5, 3.5, 3.5), nrow=n, ncol=K, byrow=T)
															
	true_result[,,4] <- matrix(c(1.7, 10.7, 10.7, 10.7,
															2.7, 2.7, 20.7, 20.7,
															3.7, 3.7, 3.7, 30.7), nrow=n, ncol=K, byrow=T)
																
	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
									 
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

	net_states <- array(NA, c(n,K,T_))
    
	net_states[,,1] <- matrix(c(0,0,0,0,
															0,0,0,0,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
    
	net_states[,,2] <- matrix(c(0,1,1,1,
                                0,0,0,0,
                                0,0,0,0), nrow=n, ncol=K, byrow=T)
	
	net_states[,,3] <- matrix(c(0,1,1,1,
															0,0,1,1,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
	
	net_states[,,4] <- matrix(c(0,1,1,1,
															0,0,1,1,
															0,0,0,1), nrow=n, ncol=K, byrow=T)
	
	
	active_mu <- matrix(c(10.1, 10.3, 10.5, 10.7,
												20.1, 20.3, 20.5, 20.7,
												30.1, 30.3, 30.5, 30.7), nrow=n, ncol=T_, byrow=T)

	active_sd <- matrix(rep(0.01, n*K), nrow=n, ncol=T_)
	
	inactive_mu <- matrix(c(1.1, 1.3, 1.5, 1.7,
													2.1, 2.3, 2.5, 2.7,
													3.1, 3.3, 3.5, 3.7), nrow=n, ncol=T_, byrow=T)

	inactive_sd <- matrix(rep(0.01, n*K), nrow=n, ncol=T_)
	
	obs_mat <- getObsMat(act_mat=NULL, net_states,  active_mu, active_sd, inactive_mu, inactive_sd, mu_type="perGeneTime")
	checkEquals(true_result, obs_mat, tolerance=(max(active_sd) + max(inactive_sd)))
}


test.getObsMatMuTypePerGeneExpTime_nodeStates <- function() {

	n <- 3
	K <- 4
	T_ <- 4
	
	true_result <- array(NA, c(n,K,T_))
    
	true_result[,,1] <- matrix(c(1.1, 1.3, 1.5, 1.7,
															 1.1, 1.3, 1.5, 1.7,
															 1.1, 1.3, 1.5, 1.7), nrow=n, ncol=K, byrow=T)
	
	true_result[,,2] <- matrix(c(2.1, 20.3, 20.5, 20.7,                                 
															2.1, 2.3, 2.5, 2.7,
															2.1, 2.3, 2.5, 2.7), nrow=n, ncol=K, byrow=T)
	
	true_result[,,3] <- matrix(c(3.1, 30.3, 30.5, 30.7,                                 
															3.1, 3.3, 30.5, 30.7,
															3.1, 3.3, 3.5, 3.7), nrow=n, ncol=K, byrow=T)
	
	true_result[,,4] <- matrix(c(4.1, 40.3, 40.5, 40.7,                                 
															4.1, 4.3, 40.5, 40.7,
															4.1, 4.3, 4.5, 40.7), nrow=n, ncol=K, byrow=T)

	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

	net_states <- array(NA, c(n,K,T_))
	
	net_states[,,1] <- matrix(c(0,0,0,0,
															0,0,0,0,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
	
	net_states[,,2] <- matrix(c(0,1,1,1,
															0,0,0,0,
															0,0,0,0), nrow=n, ncol=K, byrow=T)

	net_states[,,3] <- matrix(c(0,1,1,1,
															0,0,1,1,
															0,0,0,0), nrow=n, ncol=K, byrow=T)
	
	net_states[,,4] <- matrix(c(0,1,1,1,
															0,0,1,1,
															0,0,0,1), nrow=n, ncol=K, byrow=T)
	
	active_mu <- array(NA, c(n,K,T_))
	
	active_mu[,,1] <- matrix(c(10.1, 10.3, 10.5, 10.7,
															 10.1, 10.3, 10.5, 10.7,
															 10.1, 10.3, 10.5, 10.7), nrow=n, ncol=K, byrow=T)
	
	active_mu[,,2] <- matrix(c(20.1, 20.3, 20.5, 20.7,                                 
															20.1, 20.3, 20.5, 20.7,
															20.1, 20.3, 20.5, 20.7), nrow=n, ncol=K, byrow=T)
	
	active_mu[,,3] <- matrix(c(30.1, 30.3, 30.5, 30.7,                                 
															30.1, 30.3, 30.5, 30.7,
															30.1, 30.3, 30.5, 30.7), nrow=n, ncol=K, byrow=T)
	
	active_mu[,,4] <- matrix(c(40.1, 40.3, 40.5, 40.7,                                 
															40.1, 40.3, 40.5, 40.7,
															40.1, 40.3, 40.5, 40.7), nrow=n, ncol=K, byrow=T)
	
	active_sd <-  array(0.01, c(n,K,T_))

	inactive_mu <- array(NA, c(n,K,T_))
	inactive_mu[,,1] <- matrix(c(1.1, 1.3, 1.5, 1.7,
															 1.1, 1.3, 1.5, 1.7,
															 1.1, 1.3, 1.5, 1.7), nrow=n, ncol=K, byrow=T)

	inactive_mu[,,2] <- matrix(c(2.1, 2.3, 2.5, 2.7,
															 2.1, 2.3, 2.5, 2.7,
															 2.1, 2.3, 2.5, 2.7), nrow=n, ncol=K, byrow=T)
	
	inactive_mu[,,3] <- matrix(c(3.1, 3.3, 3.5, 3.7,
															 3.1, 3.3, 3.5, 3.7,
															 3.1, 3.3, 3.5, 3.7), nrow=n, ncol=K, byrow=T)

	inactive_mu[,,4] <- matrix(c(4.1, 4.3, 4.5, 4.7,
															 4.1, 4.3, 4.5, 4.7,
															 4.1, 4.3, 4.5, 4.7), nrow=n, ncol=K, byrow=T)

	inactive_sd <- array(0.01, c(n,K,T_))
	
	obs_mat <- getObsMat(act_mat=NULL, net_states,  active_mu, active_sd, inactive_mu, inactive_sd, mu_type="perGeneExpTime")
	checkEquals(true_result, obs_mat, tolerance=(max(active_sd) + max(inactive_sd)))
}
