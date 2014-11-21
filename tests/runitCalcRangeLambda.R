test.calcRangeLambda <- function() {

	n <- 3
	K <- 4
	
	true_result <- c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
									 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.25)

	
	obs_mat <- matrix(c(0.56, 0.95, 0.95, 0.95,
											0.56, 0.56, 0.95, 0.95,
											0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)

	delta <- rep(0.755, n)
	delta_type <- "perGene"
	
	lambda <- calcRangeLambda(obs=obs_mat, delta=delta, delta_type=delta_type)
	
	checkEquals(true_result, lambda)
}


test.calcRangeLambdaPerGeneExp<- function() {

	n <- 3
	K <- 4

	true_result <- c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 
									 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.33)

	
	obs_mat <- matrix(c(0.56, 0.95, 0.95, 0.95,
											0.56, 0.56, 0.95, 0.95,
											0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)

	delta = matrix(c(0.755, 0.755, 0.96, 0.755, 
									 0.755, 0.755, 0.96, 0.755,
									 0.755, 0.755, 0.96, 0.755), nrow=n, ncol=K, byrow=TRUE)
	delta_type <- "perGeneExp"
	
	lambda <- calcRangeLambda(obs=obs_mat, delta=delta, delta_type=delta_type)
	
	checkEquals(true_result, lambda)
}


test.calcRangeLambdaTimeSeries <- function() {

	n <- 3
	K <- 4
	T_ <- 4

	true_result <- c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
									 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 
									 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 
									 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 
									 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 
									 0.90, 0.92, 0.94, 0.96, 0.98, 1.00, 1.05, 1.09)
	
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

	delta <- rep(0.755, n)
	delta_type <- "perGene"
	
	lambda <- calcRangeLambda(obs=obs_mat, delta=delta, delta_type=delta_type, flag_time_series=TRUE)
	
	checkEquals(true_result, lambda)
}

test.calcRangeLambdaTimeSeriesPerGeneExp <- function() {

	n <- 3
	K <- 4
	T_ <- 4

	true_result <- c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
									 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 
									 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 
									 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 
									 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 
									 0.90, 0.92, 0.94, 0.96, 0.98, 1.00, 1.05, 1.10, 1.15, 1.20,
									 1.25, 1.28)
	
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

	delta = matrix(c(0.755, 0.755, 0.96, 0.755, 
									 0.755, 0.755, 0.96, 0.755,
									 0.755, 0.755, 0.96, 0.96), nrow=n, ncol=K, byrow=TRUE)
	delta_type <- "perGeneExp"
	
	lambda <-calcRangeLambda(obs=obs_mat, delta=delta, delta_type=delta_type, flag_time_series=TRUE)
	
	checkEquals(true_result, lambda)
}


test.calcRangeLambdaTimeSeriesPerGeneTime <- function() {

	n <- 3
	K <- 4
	T_ <- 4

	true_result <- c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
									 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 
									 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 
									 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 
									 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 
									 0.90, 0.92, 0.94, 0.96, 0.98, 1.00, 1.05, 1.10, 1.15, 1.20,
									 1.25)
	
	obs_mat = array(NA, c(n,K,T_))

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

	delta <- matrix(c(0.755, 0.755, 0.96, 0.755, 
									 0.755, 0.755, 0.96, 0.755,
									 0.755, 0.755, 0.96, 0.755), nrow=n, ncol=T_, byrow=TRUE)
	delta_type <- "perGeneTime"
	
	lambda <- calcRangeLambda(obs=obs_mat, delta=delta, delta_type=delta_type, flag_time_series=TRUE)
	
	checkEquals(true_result, lambda)
}


test.calcRangeLambdaTimeSeriesperGeneExpTime <- function() {

	n <- 3
	K <- 4
	T_ <- 4

	true_result <- c(0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
									 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 
									 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 
									 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 
									 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.88, 
									 0.90, 0.92, 0.94, 0.96, 0.98, 1.00, 1.05, 1.10, 1.15, 1.19)
	
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

	delta <- array(NA, c(n,K,T_))
	
	delta[,,1] <- matrix(c(0.755, 0.755, 0.96, 0.755, 
											  0.755, 0.755, 0.96, 0.755,
											  0.755, 0.755, 0.96, 0.755), nrow=n, ncol=K, byrow=TRUE)
											  
	delta[,,2] <- matrix(c(0.755, 0.755, 0.96, 0.755, 
											  0.755, 0.755, 0.96, 0.755,
											  0.755, 0.755, 0.96, 0.755), nrow=n, ncol=K, byrow=TRUE)
											  
	delta[,,3] <- matrix(c(0.755, 0.755, 0.755, 0.755, 
											  0.755, 0.755, 0.755, 0.755,
											  0.755, 0.755, 0.755, 0.755), nrow=n, ncol=K, byrow=TRUE)
											  
	delta[,,4] <- matrix(c(0.755, 0.755, 0.96, 0.755, 
											  0.755, 0.755, 0.96, 0.755,
											  0.755, 0.755, 0.96, 0.755), nrow=n, ncol=K, byrow=TRUE)
									 
	delta_type <- "perGeneExpTime"
	
	lambda <- calcRangeLambda(obs=obs_mat, delta=delta, delta_type=delta_type, flag_time_series=TRUE)
	
	checkEquals(true_result, lambda)
}
