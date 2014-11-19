#
# generate observation matrix for simulated data, either from activation matrix or network states
#
getObsMat <- function(act_mat=NULL, net_states=NULL, active_mu, active_sd, inactive_mu, inactive_sd, mu_type) {
  
	if (!is.null(net_states)) {
		T_ <- dim(net_states)[3]
		obs_mat <- net_states
		n <- dim(obs_mat)[1]
		K <- dim(obs_mat)[2]

		if (mu_type == "single") {
			for (t in 1:T_){
				obs_mat[,,t] <- .set_single_values(obs_mat[,,t], net_states[,,t], active_mu, active_sd, inactive_mu, inactive_sd)
			}
		}
		else if (mu_type == "perGene") {
			for (t in 1:T_){
				obs_mat[,,t] <- .set_per_gene_values(obs_mat[,,t], net_states[,,t], active_mu, active_sd, inactive_mu, inactive_sd, n)
			}
		}
		else if (mu_type == "perGeneExp") {
			for (t in 1:T_){
				obs_mat[,,t] <- .set_per_gene_exp_values(obs_mat[,,t], net_states[,,t], active_mu, active_sd, inactive_mu, inactive_sd, n, K)
			}
		}
		else if (mu_type == "perGeneTime") {
			for (t in 1:T_){
				obs_mat[,,t] <- .set_per_gene_time_values(obs_mat[,,t], net_states[,,t], active_mu, active_sd, inactive_mu, inactive_sd, n, t)
			}
		}
		else if (mu_type == "perGeneExpTime") {
			for (t in 1:T_){
				obs_mat[,,t] <- .set_per_gene_exp_time_values(obs_mat[,,t], net_states[,,t], active_mu, active_sd, inactive_mu, inactive_sd, n, K, t)
			}
		}
	}
	else {
        obs_mat <- act_mat
        n <- dim(obs_mat)[1]
				K <- dim(obs_mat)[2]
        if (mu_type == "single") {
            obs_mat <- .set_single_values(obs_mat, act_mat, active_mu, active_sd, inactive_mu, inactive_sd)
        }
        else if (mu_type == "perGene") {    
            obs_mat <- .set_per_gene_values(obs_mat, act_mat, active_mu, active_sd, inactive_mu, inactive_sd, n)
        }
        else if (mu_type == "perGeneExp") {
            obs_mat <- .set_per_gene_exp_values(obs_mat, act_mat, active_mu, active_sd, inactive_mu, inactive_sd, n, K)
        }
    }
  
  return(obs_mat)
}


.set_single_values <- function(obs_mat, act_mat, active_mu, active_sd, inactive_mu, inactive_sd) {

	obs_mat[act_mat == 1] <- rnorm(sum(act_mat == 1), active_mu, active_sd)
	obs_mat[act_mat == 0] <- rnorm(sum(act_mat == 0), inactive_mu, inactive_sd)
	
	return(obs_mat)
}


.set_per_gene_values <- function(obs_mat, act_mat, active_mu, active_sd, inactive_mu, inactive_sd, n) {

	for (i in 1:n) {
		obs_mat[i, ][act_mat[i, ] == 1] <- rnorm(sum(act_mat[i, ] == 1), active_mu[i], active_sd[i])
		obs_mat[i, ][act_mat[i, ] == 0] <- rnorm(sum(act_mat[i, ] == 0), inactive_mu[i], inactive_sd[i])
	}
	
	return(obs_mat)
}


.set_per_gene_exp_values <- function(obs_mat, act_mat, active_mu, active_sd, inactive_mu, inactive_sd, n, K) {

	for (i in 1:n) {
		for (k in 1:K) {
			if (act_mat[i,k] == 1) {
				obs_mat[i,k] <- rnorm(1, active_mu[i,k], active_sd[i,k])
			}
			else {
				obs_mat[i,k] <- rnorm(1, inactive_mu[i,k], inactive_sd[i,k])
			}
		}
	}
	
	return(obs_mat)
}


.set_per_gene_time_values <- function(obs_mat, act_mat, active_mu, active_sd, inactive_mu, inactive_sd, n, t) {

	for (i in 1:n) {
		obs_mat[i, ][act_mat[i, ] == 1] <- rnorm(sum(act_mat[i, ] == 1), active_mu[i,t], active_sd[i,t])
		obs_mat[i, ][act_mat[i, ] == 0] <- rnorm(sum(act_mat[i, ] == 0), inactive_mu[i,t], inactive_sd[i,t])
    }
	
	return(obs_mat)
}


.set_per_gene_exp_time_values <- function(obs_mat, act_mat, active_mu, active_sd, inactive_mu, inactive_sd, n, K, t) {

	for (i in 1:n) {
		for (k in 1:K) {
			if (act_mat[i,k] == 1) {
				obs_mat[i,k] <- rnorm(1, active_mu[i,k,t], active_sd[i,k,t])
			}
			else {
				obs_mat[i,k] <- rnorm(1, inactive_mu[i,k,t], inactive_sd[i,k,t])
			}
		}
	}
	
	return(obs_mat)
}
