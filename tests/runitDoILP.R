.setUp <- function(){
	
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
											
	lambda <<- 1/10
	annot <<- getEdgeAnnot(n)
}


test.doILPShortExamplePerGene <- function() {

	true_result_objval <- 13.52785
	true_result_solution <- c(0.0000000, 0.7947368, 0.0000000, 
														0.0000000, 0.0000000, 1.9358974,
														0.0000000, 0.0000000, 0.0000000,
														0.0000000, 0.0000000, 1.1411606,
														0.0000000, 0.0000000, 0.0000000,
														0.0000000, 0.0000000, 0.0000000,
														0.7550000, 0.0000000, 0.0000000,
														0.0000000, 0.4450526, 0.4450526,
														0.0000000, 0.0000000, 0.0000000,
														0.0000000, 0.0000000, 0.0000000,
														0.0000000, 0.0000000, 0.0000000)
	
	delta = rep(0.755, n)
	delta_type <- "perGene"
	
	res <- doILP(obs_mat, delta, lambda, b, n, K, T_=NULL, annot, delta_type, prior=NULL, sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE)

	checkEquals(true_result_objval, res$objval, tolerance=0.00001)
	checkEquals(true_result_solution, res$solution, tolerance=0.00001)
}


test.doILPShortExamplePerGeneExp <- function() {

	true_result_objval <- 19.68196
	true_result_solution <- c(0.0000000, 0.0000000, 0.0000000,
													  0.0000000, 0.0000000, 1.9358974,
													  1.9358974, 1.9358974, 0.0000000,
													  0.0000000, 1.1411606, 1.1411606,
													  1.9358974, 0.0000000, 0.0000000,
													  0.0000000, 0.0000000, 0.0000000,
													  0.7550000, 0.0000000, 0.0000000,
														0.0000000, 0.4450526, 0.4450526,
														0.0000000, 0.0000000, 0.0000000,
														0.0000000, 0.0000000, 0.0000000,
														0.0000000, 0.0000000, 0.0000000)

	delta = matrix(c(0.755, 0.755, 0.96, 0.755, 
									 0.755, 0.755, 0.96, 0.755,
									 0.755, 0.755, 0.96, 0.755), nrow=n, ncol=K, byrow=TRUE)

	delta_type <- "perGeneExp"
	
	res <- doILP(obs_mat, delta, lambda, b, n, K, T_=NULL, annot, delta_type, prior=NULL, sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE)

	checkEquals(true_result_objval, res$objval, tolerance=0.00001)
	checkEquals(true_result_solution, res$solution, tolerance=0.00001)
}

