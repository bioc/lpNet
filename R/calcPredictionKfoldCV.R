#
# predict node states + observation values for KfoldCV
#
calcPredictionKfoldCV <- function(obs, delta, b, n, K, adja, baseline, rem_entries=NULL, rem_entries_vec=NULL,
																 active_mu, active_sd, inactive_mu, inactive_sd, mu_type, flag_time_series=FALSE) {
 
  if (flag_time_series == FALSE) {
		predict <- .calcPredictionKfoldCV_steadyState(obs, delta, b, n, K, adja, baseline,
																								  active_mu, active_sd, inactive_mu, inactive_sd, mu_type)
		
  }
  else {
		predict <- .calcPredictionKfoldCV_timeSeries(obs, delta, b ,n ,K, adja, baseline, rem_entries, rem_entries_vec,
																								 active_mu, active_sd, inactive_mu, inactive_sd, mu_type)
  }
  
  return(predict)
}


#
# kfoldCV prediction for LP original
#
.calcPredictionKfoldCV_steadyState <- function(obs, delta, b, n, K, adja, baseline,
																							 active_mu, active_sd, inactive_mu, inactive_sd, mu_type) {

  act_mat <- calcActivation(adja, b, n, K)
  predict <- getObsMat(act_mat, net_states=NULL, active_mu, active_sd, inactive_mu, inactive_sd, mu_type)
  
  return(predict)
}


#
# kfoldCV prediction for half discretized model
#
.calcPredictionKfoldCV_timeSeries <- function(obs, delta, b ,n ,K, adja, baseline, rem_entries, rem_entries_vec,
																						  active_mu, active_sd, inactive_mu, inactive_sd, mu_type) {

	inact_entries <- which(b == 0)
	predict <- obs

	if (dim(rem_entries)[1] > 0) {
		if (mu_type == "single") {
		
			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <- rem_entries[ent,3]
				
				rem_ent_test <- rem_entries_vec[ent] %% (n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K  # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- rem_ent_test %in% inact_entries
				
				delta_rem <- delta[rem_gene]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b, active_mu, active_sd, 
																											inactive_mu, inactive_sd, rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		else if (mu_type == "perGene") {

			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <- rem_entries[ent,3]
				
				rem_ent_test <- rem_entries_vec[ent] %% (n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K  # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- rem_ent_test %in% inact_entries
				
				delta_rem <- delta[rem_gene]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b, active_mu[rem_gene], 
																											active_sd[rem_gene], inactive_mu[rem_gene], inactive_sd[rem_gene], 
																											rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		else if (mu_type == "perGeneExp") {
		
			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <-rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent] %% (n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- rem_ent_test %in% inact_entries
				
				delta_rem <- delta[rem_gene, rem_k]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b, active_mu[rem_gene, rem_k], 
																											active_sd[rem_gene, rem_k], inactive_mu[rem_gene, rem_k], inactive_sd[rem_gene, rem_k], 
																											rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		
		else if (mu_type == "perGeneTime") {
		
			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <- rem_entries[ent,3]
				
				rem_ent_test <- rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K  # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- (rem_ent_test %in% inact_entries)
				
				delta_rem <- delta[rem_gene, rem_t]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b,
																											active_mu[rem_gene, rem_t], active_sd[rem_gene, rem_t], 
																											inactive_mu[rem_gene, rem_t], inactive_sd[rem_gene, rem_t], 
																											rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		
		else if (mu_type == "perGeneExpTime") {
		
			for (ent in 1:dim(rem_entries)[1]) {
				rem_gene <- rem_entries[ent,1]
				rem_k <- rem_entries[ent,2]
				rem_t <- rem_entries[ent,3]
				
				rem_ent_test = rem_entries_vec[ent]%%(n*K)
				if (rem_ent_test == 0) 
					rem_ent_test <- n*K  # when the removed entry is in the last position the remainder of the division is zero and the test fails
				res <- rem_ent_test %in% inact_entries
				
				delta_rem <- delta[rem_gene, rem_k, rem_t]
				predict <- .calculatePredictionValue_Kfold_ts(predict, obs, n, adja, baseline, delta, delta_rem, b,
																											active_mu[rem_gene, rem_k, rem_t], active_sd[rem_gene, rem_k, rem_t], 
																											inactive_mu[rem_gene, rem_k, rem_t], inactive_sd[rem_gene, rem_k, rem_t], 
																											rem_gene, rem_k, rem_t, res, mu_type)
			}
		}
		
	}
	return(predict)
}


.calculatePredictionValue_Kfold_ts <- function(predict, obs, n, adja, baseline, delta, delta_rem, b, active_mu, active_sd, 
																							 inactive_mu, inactive_sd, rem_gene, rem_k, rem_t, res, mu_type) {

	# if the removed entry is an inactive node due to some knockdown, then predict as inactivet=rem_entries[ent,3]
	if (res == TRUE) {
		predict[rem_gene,rem_k,rem_t] <- rnorm(1, inactive_mu, inactive_sd)
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
				predict[rem_gene,rem_k,rem_t] <- rnorm(1, active_mu, active_sd) 
			}
			else {
				predict[rem_gene,rem_k,rem_t] <- rnorm(1, inactive_mu, inactive_sd)
			}
		}
		else {
			flagNA <- 0
			for (j in 1:length(pa)) {
				
				delta_pa <- .setDeltaValue_KfoldCV_ts(delta, pa, j, rem_k, rem_t, mu_type)
				if (is.na(obs[pa[j],rem_k,rem_t-1])) {  # if parent observation is NA
					flagNA <- 1
					
					if ((adja[pa[j],rem_gene] < 0) & (b[(rem_k-1)*n + pa[j]] == 1)) {  # if the incoming edge is negative and the parent node is active, node state is unknown
						predict[rem_gene,rem_k,rem_t] <- NA
						return(predict=predict)
					}
				}
				else if ((obs[pa[j],rem_k,rem_t-1] >= delta_pa) & (b[(rem_k-1)*n + pa[j]] == 1)) {  # if parent is active and not silenced calculate node inflow
					in_flow <- sum(in_flow, adja[pa[j],rem_gene] * obs[pa[j],rem_k,rem_t-1], na.rm=T)
				}
			}
			
			if ((flagNA == 1) & (in_flow < delta_rem)) {  # if parent observation is NA and inflow < delta, node state is unknown
				predict[rem_gene,rem_k,rem_t] <- NA
				return(predict)
			}			
			
			# predict node state according to inflow
			if (in_flow >= delta_rem) {  
				predict[rem_gene,rem_k,rem_t] <- rnorm(1, active_mu, active_sd) 
			}
			else {
				predict[rem_gene,rem_k,rem_t] <- rnorm(1, inactive_mu, inactive_sd)
			}
		}
	}
	return(predict)
}


.setDeltaValue_KfoldCV_ts <- function(delta, pa, j, rem_k, rem_t, mu_type) {

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
