#
# predict node state + observation value for LOOCV
#
calcPredictionLOOCV <-function(obs, delta, b, n ,K, adja, baseline, rem_gene, rem_k, rem_t=NULL, active_mu, 
															 active_sd, inactive_mu, inactive_sd, mu_type, flag_time_series=FALSE) {
	
	if (flag_time_series == FALSE){
		predict <- .calcPredictionLOOCV_steadyState(obs, delta, b, n ,K, adja, baseline, rem_gene, rem_k, 
																							  active_mu, active_sd, inactive_mu, inactive_sd, mu_type)
	}
	else {
		predict <- .calcPredictionLOOCV_timeSeries(obs, delta, b, n, K, adja, baseline, rem_gene, rem_k, rem_t,
																							 active_mu, active_sd, inactive_mu, inactive_sd, mu_type)
	}
	
	return(predict)
}


.calcPredictionLOOCV_steadyState <- function(obs, delta, b, n ,K, adja, baseline, rem_gene, rem_k, 
																						 active_mu, active_sd, inactive_mu, inactive_sd, mu_type) {
	
	kds <- matrix(b, nrow=n, ncol=K)
	sil_gene_ids <- which(kds[ ,rem_k] == 0)

	if (mu_type == "single") {
		predict <- .calculatePredictionValue_LOOCV_ss(predict, obs, delta[rem_gene], adja, active_mu, 
																									active_sd, inactive_mu,
																									inactive_sd,rem_gene, rem_k, sil_gene_ids)
	}
	else if (mu_type == "perGene") {
		predict <- .calculatePredictionValue_LOOCV_ss(predict, obs, delta[rem_gene], adja, active_mu[rem_gene], 
																									active_sd[rem_gene], inactive_mu[rem_gene],
																									inactive_sd[rem_gene],rem_gene, rem_k, sil_gene_ids)
	}
	else if (mu_type == "perGeneExp"){
		predict <- .calculatePredictionValue_LOOCV_ss(predict, obs, delta[rem_gene, rem_k], adja, active_mu[rem_gene, rem_k], 
																									active_sd[rem_gene, rem_k], inactive_mu[rem_gene, rem_k],
																									inactive_sd[rem_gene, rem_k],rem_gene, rem_k, sil_gene_ids)
	}
	return(predict)
}

.calculatePredictionValue_LOOCV_ss <- function(predict, obs, delta, adja, active_mu, active_sd, 
																							 inactive_mu, inactive_sd, rem_gene, rem_k, sil_gene_ids) {
																							
	if (rem_gene %in% sil_gene_ids) {
		predict <- rnorm(1, inactive_mu, inactive_sd)
	}
	else {  # else: in_flow: sum of all parents of rem_gene after knockdown rem_k times the weights
		pa <- which(adja[ ,rem_gene] != 0)
					
		if(length(pa)==0){  # if there are no parents: rem_gene is root node and thus is considered to be active since it is not silenced
			predict <- rnorm(1,active_mu,active_sd)
		}
		else { # calculate in_flow
			in_flow <- 0
			for (j in 1:length(pa)) {
				in_flow <- sum(in_flow, adja[pa[j],rem_gene] * obs[pa[j],rem_k], na.rm=T)
			}
			if (in_flow >= delta) {
				predict <- rnorm(1, active_mu, active_sd)
			}
			else {
				predict <- rnorm(1, inactive_mu, inactive_sd)
			}
		}
	}
	
	return(predict)
}


.calcPredictionLOOCV_timeSeries <- function(obs, delta, b, n, K, adja, baseline, rem_gene, rem_k, rem_t,
																						active_mu, active_sd, inactive_mu, inactive_sd, mu_type) { 
	
	kds <- matrix(b, nrow=n, ncol=K)
	sil_gene_ids <- which(kds[ ,rem_k] == 0)
	
	if (mu_type == "single") {
		delta_rem <- delta[rem_gene]
		predict <- .calculatePredictionValue_LOOCV_ts(predict, obs, delta, delta_rem, b, n, adja, baseline, active_mu, active_sd, 
																									inactive_mu, inactive_sd, rem_gene, rem_k, rem_t, sil_gene_ids, mu_type)
	}
	# if there is an in/active_mu/sd per gene
	else if (mu_type == "perGene") {
		delta_rem <- delta[rem_gene]
		predict <- .calculatePredictionValue_LOOCV_ts(predict, obs, delta, delta_rem, b, n, adja, baseline, active_mu[rem_gene], active_sd[rem_gene], 
																									inactive_mu[rem_gene], inactive_sd[rem_gene], rem_gene, rem_k, rem_t, sil_gene_ids, mu_type)
	}
	# if there is an in/active_mu/sd and delta per gene per knockdown exp
	else if (mu_type == "perGeneExp") {
		delta_rem <- delta[rem_gene, rem_k]
		predict <- .calculatePredictionValue_LOOCV_ts(predict, obs, delta, delta_rem, b, n, adja, baseline, active_mu[rem_gene, rem_k], 
																									active_sd[rem_gene, rem_k], inactive_mu[rem_gene, rem_k], inactive_sd[rem_gene, rem_k],
																									rem_gene, rem_k, rem_t, sil_gene_ids, mu_type)
	}
	# if there is an in/active_mu/sd and delta per gene per time point
	else if (mu_type == "perGeneTime") {
		delta_rem <- delta[rem_gene, rem_t]
		predict <- .calculatePredictionValue_LOOCV_ts(predict, obs, delta, delta_rem, b, n, adja, baseline, active_mu[rem_gene, rem_t],
																									active_sd[rem_gene, rem_t], inactive_mu[rem_gene, rem_t], inactive_sd[rem_gene, rem_t],
																									rem_gene, rem_k, rem_t, sil_gene_ids, mu_type)
	}
	else if (mu_type == "perGeneExpTime") {
		delta_rem <- delta[rem_gene, rem_k, rem_t]
		predict <- .calculatePredictionValue_LOOCV_ts(predict, obs, delta, delta_rem, b, n, adja, baseline, active_mu[rem_gene, rem_k, rem_t],
																									active_sd[rem_gene, rem_k, rem_t], inactive_mu[rem_gene, rem_k, rem_t], 
																									inactive_sd[rem_gene, rem_k, rem_t], rem_gene, rem_k, rem_t, sil_gene_ids, mu_type)
	}
	
	
	return(predict)
}
	

.calculatePredictionValue_LOOCV_ts <- function(predict, obs, delta, delta_rem, b, n, adja, baseline, active_mu, active_sd, 
																							inactive_mu, inactive_sd, rem_gene, rem_k, rem_t, sil_gene_ids, mu_type) {

	if (rem_gene %in% sil_gene_ids) {  # if the removed entry is an inactive node due to some knockdown, then predict as inactive
		predict <- rnorm(1, inactive_mu, inactive_sd)
	}
	else {
		if (is.na(baseline[rem_gene])) {
			in_flow <- 0
		}
		else {
			in_flow <- baseline[rem_gene]
		}
		
		pa <- which(adja[ ,rem_gene] != 0)
		if (length(pa) == 0) {  # if there are no parents: rem_gene is root node 
			if (in_flow >= delta_rem) {  # root node is active if its inflow is greater than its delta
				predict <- rnorm(1, active_mu, active_sd) 
			}
			else {
				predict <- rnorm(1, inactive_mu, inactive_sd)
			}
		}
		else {
			flagNA <- 0
			for (j in 1:length(pa)) {

				delta_pa <- .setDeltaValue_LOOCV_ts(delta, pa, j, rem_k, rem_t, mu_type)
				if (is.na(obs[pa[j],rem_k,rem_t-1])) {  # if parent observation is NA
					flagNA <- 1
					
					if ((adja[pa[j],rem_gene] < 0) & (b[(rem_k-1)*n + pa[j]] == 1)) {  # if the incoming edge is negative and the parent node is active, node state is unknown
						predict <- NA
						return(predict)
					}
				}
				else if ((obs[pa[j],rem_k,rem_t-1] >= delta_pa) & (b[(rem_k-1)*n + pa[j]] == 1)) { # if parent is active and not silenced calculate node inflow
					in_flow <- sum(in_flow, adja[pa[j],rem_gene] * obs[pa[j],rem_k,rem_t-1], na.rm=T)
				}
			}
			
			if ((flagNA == 1) & (in_flow < delta_rem)) {  # if parent observation is NA and inflow < delta, node state is unknown
				predict <- NA
				return(predict)
			}			
			
			# predict node state according to inflow
			if (in_flow >= delta_rem) {
				predict<- rnorm(1, active_mu, active_sd) 
			}
			else {
				predict <- rnorm(1, inactive_mu, inactive_sd)
			}
		}
	}
	return(predict)
}


.setDeltaValue_LOOCV_ts <- function(delta, pa, j, rem_k, rem_t, mu_type) {

	if ((mu_type == "single") | (mu_type == "perGene")) {
		delta_pa <- delta[pa[j]]
	}
	else if (mu_type == "perGeneExp") {
		delta_pa <- delta[pa[j], rem_k]
	}
	else if (mu_type == "perGeneTime") {
		delta_pa <- delta[pa[j], rem_t-1]
	}
	else if (mu_type == "perGeneExpTime") {
		delta_pa <- delta[pa[j], rem_k, rem_t-1]
	}

	return(delta_pa)
}
