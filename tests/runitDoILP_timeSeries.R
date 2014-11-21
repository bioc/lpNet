.setUp <- function() {

	n <<- 3
	K <<- 4
	T_ <<- 4
	
	T_nw <<- matrix(c(0,1,0,
									 0,0,1,
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
														0.56, 0.56, 0.56, 0.56,
														0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)

	obs_mat[,,3] <<- matrix(c(0.56, 0.95, 0.95, 0.95,
														0.56, 0.56, 0.95, 0.95,
														0.56, 0.56, 0.56, 0.56), nrow=n, ncol=K, byrow=TRUE)
	
	obs_mat[,,4] <<- matrix(c(0.56, 0.95, 0.95, 0.95,
														0.56, 0.56, 0.95, 0.95,
														0.56, 0.56, 0.56, 0.95), nrow=n, ncol=K, byrow=TRUE)
														
	lambda <<- 1/10
	annot <<- getEdgeAnnot(n)
}


test.doILPTimeSeriesShortExamplePerGene <- function() {
	
	true_result_objval <- 2.344474
	true_result_solution <- c(0.0000000, 0.7947368, 0.0000000, 
													  0.0000000, 0.0000000, 0.7947368, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.7550000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000,
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000)
	
	delta <- rep(0.755, n)

	delta_type <- "perGene"

	res <- doILP(obs_mat, delta, lambda, b, n, K, T_, annot, delta_type, prior=NULL, 
							 sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE, flag_time_series=TRUE)
	
	checkEquals(true_result_objval, res$objval, tolerance=0.00001)
	checkEquals(true_result_solution, res$solution, tolerance=0.00001)
}


test.doILPTimeSeriesShortExamplePerGenePerExp <- function() {


	true_result_objval <- 24.99447
	true_result_solution <- c(0.0000000, 0.7947368, 0.0000000, 
													  0.0000000, 0.0000000, 0.7947368, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.7550000, 0.0000000, 0.0000000,
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.7550000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.7550000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000, 
													  0.7550000, 0.0000000, 0.0000000, 
													  0.0000000, 0.0000000, 0.0000000)

	delta <- matrix(c(0.755, 0.755, 0.96, 0.755, 
									  0.755, 0.755, 0.96, 0.755,
									  0.755, 0.755, 0.96, 0.755), nrow=n, ncol=K, byrow=TRUE)
									 
	delta_type <- "perGeneExp"
	
	res <- doILP(obs_mat, delta, lambda, b, n, K, T_, annot, delta_type, prior=NULL, 
							sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE, flag_time_series=TRUE)
		
	checkEquals(true_result_objval, res$objval, tolerance=0.00001)
	checkEquals(true_result_solution, res$solution, tolerance=0.00001)
}


test.doILPTimeSeriesShortExamplePerGenePerTime <- function() {


	true_result_objval <- 109.5545
	true_result_solution <- c(0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.7947368, 0.7947368, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.7550000, 0.7550000, 0.7550000, 
														0.0000000, 0.7550000, 0.7550000, 
														0.0000000, 0.0000000, 0.7550000, 
														0.0000000, 0.7550000, 0.0000000, 
														0.0000000, 0.7550000, 0.7550000, 
														0.0000000, 0.7550000, 0.7550000, 
														0.7550000, 0.0000000, 0.0000000, 
														0.7550000, 0.0000000, 0.0000000, 
														0.7550000, 0.0000000, 0.0000000, 
														0.0000000, 0.7550000, 0.7550000, 
														0.0000000, 0.0000000, 0.7550000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000)
	
	delta <- matrix(c(0.755, 0.755, 0.96, 0.755, 
									  0.755, 0.755, 0.96, 0.755,
									  0.755, 0.755, 0.96, 0.755), nrow=n, ncol=K, byrow=TRUE)
									 
	delta_type <- "perGeneTime"
	
	res <- doILP(obs_mat, delta, lambda, b, n, K, T_, annot, delta_type, prior=NULL, 
							 sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE, flag_time_series=TRUE)
	
	checkEquals(true_result_objval, res$objval, tolerance=0.00001)
	checkEquals(true_result_solution, res$solution, tolerance=0.00001)
}

test.doILPTimeSeriesShortExamplePerGenePerExpPerTime <- function() {

	true_result_objval <- 62.70474
	true_result_solution <- c(0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.7947368, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.7550000, 0.7550000, 0.0000000, 
														0.0000000, 0.7550000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.7550000, 0.7550000, 0.0000000, 
														0.0000000, 0.7550000, 0.0000000, 
														0.0000000, 0.7550000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.0000000, 0.7550000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000, 
														0.7550000, 0.7550000, 0.0000000, 
														0.0000000, 0.0000000, 0.0000000)

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
	
	res <- doILP(obs_mat, delta, lambda, b, n, K, T_, annot, delta_type, prior=NULL, 
							 sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE, flag_time_series=TRUE)
													
	checkEquals(true_result_objval, res$objval, tolerance=0.00001)
	checkEquals(true_result_solution, res$solution, tolerance=0.00001)
}

