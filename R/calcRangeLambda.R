#
# calculate possible lambda values. 
#   idea: the possible values increase exponentially
#
calcRangeLambda <- function(obs, delta, delta_type, flag_time_series=FALSE) {
 
	if (flag_time_series == FALSE){
		lambda_values <- .calcRangeLambda_steadyState(obs, delta, delta_type)
	}
	else {
		lambda_values <- .calcRangeLambda_timeSeries(obs,delta, delta_type)
	}
	
	return(lambda_values)
}


.calcRangeLambda_steadyState <- function(obs, delta, delta_type) {

  n_slacks <- 0
 
  if (delta_type == "perGene") {
		for (i in 1:length(delta))
			n_slacks <- sum(c(n_slacks, obs[i,] < delta[i]), na.rm=T)
	}
	else if (delta_type == "perGeneExp") {
		for (i in 1:dim(delta)[1]) {
			for (k in 1:dim(delta)[2])
				n_slacks <- sum(c(n_slacks, obs[i,k] < delta[i,k]), na.rm=T)
		}
	}
	
  if (n_slacks == 0) 
		n_slacks <- 1
  max_value <- round(n_slacks * var(c(obs), na.rm=T), digits=2)
  
	lambda <- c()
	
	if (max_value >= 0.1) {
		lambda <- c(seq(0, 0.1, by=0.01))
	}
	else {
		lambda <- c(seq(0, max_value, by=0.01))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value >= 1) {
		lambda <- c(lambda, seq(0.1, 1, by=0.02))
	}
	else {
		lambda <- c(lambda, seq(0.1, max_value, by=0.02))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value >= 2) {
		lambda <- c(lambda, seq(1, 2, by=0.05))
	}
	else {
		lambda <- c(lambda, seq(1, max_value, by=0.05))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value >= 10) {
		lambda <- c(lambda, seq(2, 10, by=1))
	}
	else{
		lambda <- c(lambda, seq(2, max_value, by=1))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value  >= 100) {
		lambda <- c(lambda, seq(10, 100, by=10))
	}
	else{
		lambda <- c(lambda, seq(10, max_value, by=10))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value  >= 1000) {
		lambda <- c(lambda, seq(100, 1000, by=100))
	}
	else{
		lambda <- c(lambda, seq(100, max_value, by=100))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value  >= 10000) {
		lambda <- c(lambda, seq(1000, 10000, by=1000))
	}
	else{
		lambda <- c(lambda, seq(1000, max_value, by=1000))
		return(unique(c(lambda, max_value)))
	}

	while (lambda[length(lambda)] < max_value){
		val = lambda[length(lambda)] * 2
		lambda <- c(lambda, val)
	}

	return(unique(c(lambda, max_value)))
}


.calcRangeLambda_timeSeries <- function(obs, delta, delta_type) {

  n_slacks <- 0
  
  if (delta_type == "perGene") {
		for (i in 1:length(delta))
			n_slacks <- sum(c(n_slacks, obs[i, , ] < delta[i]), na.rm=T)
	}
	else if (delta_type == "perGeneExp") {
		for (i in 1:dim(delta)[1]) {
			for (k in 1:dim(delta)[2])
				n_slacks <- sum(c(n_slacks, obs[i,k, ] < delta[i,k]), na.rm=T)
		}
	}
	else if (delta_type == "perGeneTime") {
		for (i in 1:dim(delta)[1]){
			for (t in 1:dim(delta)[2])
				n_slacks <- sum(c(n_slacks, obs[i, ,t] < delta[i,t]), na.rm=T)
		}
	}
	else if (delta_type == "perGeneExpTime") {
		for (i in 1:dim(delta)[1]) {
			for (k in 1:dim(delta)[2]) {
				for (t in 1:dim(delta)[3])
					n_slacks <- sum(c(n_slacks, obs[i,k,t] < delta[i,k,t]), na.rm=T)
			}
		}
	}
	
  if( n_slacks==0) 
		n_slacks <- 1
  max_value <- round(n_slacks * var(c(obs), na.rm=T), digits=2)
	
	lambda <- c()
	
	if (max_value >= 0.1) {
		lambda <- c(seq(0, 0.1, by=0.01))
	}
	else {
		lambda <- c(seq(0, max_value, by=0.01))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value >= 1) {
		lambda <- c(lambda, seq(0.1, 1, by=0.02))
	}
	else {
		lambda <- c(lambda, seq(0.1, max_value, by=0.02))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value >= 2) {
		lambda <- c(lambda, seq(1, 2, by=0.05))
	}
	else {
		lambda <- c(lambda, seq(1, max_value, by=0.05))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value >= 10) {
		lambda <- c(lambda, seq(2, 10, by=1))
	}
	else {
		lambda <- c(lambda, seq(2, max_value, by=1))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value >= 100) {
		lambda <- c(lambda, seq(10, 100, by=10))
	}
	else {
		lambda <- c(lambda, seq(10, max_value, by=10))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value >= 1000) {
		lambda <- c(lambda, seq(100, 1000, by=100))
	}
	else{
		lambda <- c(lambda, seq(100, max_value, by=100))
		return(unique(c(lambda, max_value)))
	}
	
	if (max_value >= 10000) {
		lambda <- c(lambda, seq(1000, 10000, by=1000))
	}
	else {
		lambda <- c(lambda, seq(1000, max_value, by=1000))
		return(unique(c(lambda, max_value)))
	}

	while (lambda[length(lambda)] < max_value) {
		val = lambda[length(lambda)] * 2
		lambda <- c(lambda, val)
	}
		
  return(unique(c(lambda, max_value)))
}
