test.calcActivationShortExample <- function() {
	n <- 3
	K <- 4
	
	true_result <- matrix(c(0,0,0,
													1,0,0,
													1,1,0,
													1,1,1), nrow=n, ncol=K)
	
	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

	act_mat <- calcActivation(T_nw, b, n, K)

	checkEquals(true_result, act_mat)
}


test.calcActivationShortExampleTimeSeries <- function() {
	n <- 3
	K <- 4
	
	true_result <- matrix(c(0,0,0,
													1,0,0,
													1,1,0,
													1,1,1), nrow=n, ncol=K)
	
	T_nw <- matrix(c(0,1,0,
									 0,0,1,
									 0,0,0), nrow=n, ncol=n, byrow=TRUE)
	b <- c(0,1,1,
				 1,0,1,
				 1,1,0,
				 1,1,1)

	act_mat <- calcActivation(T_nw, b, n, K, flag_gen_data=TRUE)

	checkEquals(true_result, act_mat)
}


test.calcActivation <- function() {
	n <- 5
	K <- 6
	
	true_result <- matrix(c(0,0,0,0,0,
													1,0,1,1,1,
													1,1,0,0,0,
													1,1,1,0,0,
													1,1,1,0,0,
													1,1,1,0,0), nrow=n, ncol=K)
	
	T_nw <- matrix(c(0,1,1,0,0,
									 0,0,0,-1,0,
									 0,0,0,1,0,
									 0,0,0,0,1,
									 0,0,0,0,0), nrow=n, ncol=n, byrow=TRUE)
									 
	b <- c(0,1,1,1,1,
				 1,0,1,1,1,
				 1,1,0,1,1,
				 1,1,1,0,1,
				 1,1,1,1,0,
				 1,1,1,1,1)

	act_mat <- calcActivation(T_nw, b, n, K)

	checkEquals(true_result, act_mat)
}


test.calcActivationTimeSeries <- function() {
	n <- 5
	K <- 6
	
	true_result <- matrix(c(0,0,0,0,0,
													1,0,1,1,1,
													1,1,0,1,1,
													1,1,1,0,0,
													1,1,1,1,0,
													1,1,1,1,1), nrow=n, ncol=K)
	
	T_nw <- matrix(c(0,1,1,0,0,
									 0,0,0,-1,0,
									 0,0,0,1,0,
									 0,0,0,0,1,
									 0,0,0,0,0), nrow=n, ncol=n, byrow=TRUE)
	b <- c(0,1,1,1,1,
				 1,0,1,1,1,
				 1,1,0,1,1,
				 1,1,1,0,1,
				 1,1,1,1,0,
				 1,1,1,1,1)

	act_mat <- calcActivation(T_nw, b, n, K, flag_gen_data=TRUE)

	checkEquals(true_result, act_mat)
}


test.calcActivationLargeExample <- function() {
	n <- 10
	K <- 11
	
	true_result <- matrix(c(0,0,0,1,1,1,1,1,1,1,
													1,0,0,1,1,1,1,1,1,1,
													1,0,0,1,1,1,1,1,1,1,
													1,1,1,0,0,0,0,0,0,0,
													1,1,1,1,0,0,0,0,0,0,
													1,1,1,1,1,0,0,0,0,0,
													1,0,0,1,1,1,0,0,0,0,
													1,0,0,1,1,1,1,0,0,0,
													1,0,0,1,1,1,1,1,0,0,
													1,0,0,1,1,1,1,1,1,0,
													1,0,0,1,1,1,1,1,1,1), nrow=n, ncol=K)
	
	T_nw <- matrix(c(0,1,0,0,0,0,0,0,0,0,
									 0,0,1,0,0,0,0,0,0,0,
									 0,0,0,0,0,0,0,0,0,0,
									 0,0,0,0,1,0,0,0,0,0,
									 0,0,0,0,0,1,0,0,0,0,
									 0,-1,0,0,0,0,1,0,0,0,
									 0,0,0,0,0,0,0,1,0,0,
									 0,0,0,0,0,0,0,0,1,0,
									 0,0,0,0,0,0,1,0,0,1,
									 0,0,0,0,0,0,0,0,0,0), nrow=n, ncol=n, byrow=TRUE)
	
	b <- c(0,1,1,1,1,1,1,1,1,1,
				 1,0,1,1,1,1,1,1,1,1,
				 1,1,0,1,1,1,1,1,1,1,
				 1,1,1,0,1,1,1,1,1,1,
				 1,1,1,1,0,1,1,1,1,1,
				 1,1,1,1,1,0,1,1,1,1,
				 1,1,1,1,1,1,0,1,1,1,
				 1,1,1,1,1,1,1,0,1,1,
				 1,1,1,1,1,1,1,1,0,1,
				 1,1,1,1,1,1,1,1,1,0,
				 1,1,1,1,1,1,1,1,1,1)

	act_mat <- calcActivation(T_nw, b, n, K)

	checkEquals(true_result, act_mat)
}


test.calcActivationLargeExampleTimeSeries <- function() {
	n <- 10
	K <- 11
	
	true_result <- matrix(c(0,1,1,1,1,1,1,1,1,1,
													1,0,0,1,1,1,1,1,1,1,
													1,1,0,1,1,1,1,1,1,1,
													1,1,1,0,0,0,0,0,0,0,
													1,1,1,1,0,0,0,0,0,0,
													1,1,1,1,1,0,0,0,0,0,
													1,1,1,1,1,1,0,0,0,0,
													1,1,1,1,1,1,1,0,0,0,
													1,1,1,1,1,1,1,1,0,0,
													1,1,1,1,1,1,1,1,1,0,
													1,1,1,1,1,1,1,1,1,1), nrow = n, ncol=K)
	
	T_nw <- matrix(c(0,1,0,0,0,0,0,0,0,0,
									0,0,1,0,0,0,0,0,0,0,
									0,0,0,0,0,0,0,0,0,0,
									0,0,0,0,1,0,0,0,0,0,
									0,0,0,0,0,1,0,0,0,0,
									0,-1,0,0,0,0,1,0,0,0,
									0,0,0,0,0,0,0,1,0,0,
									0,0,0,0,0,0,0,0,1,0,
									0,0,0,0,0,0,1,0,0,1,
									0,0,0,0,0,0,0,0,0,0), nrow=n, ncol=n, byrow=TRUE)
	
	b <- c(0,1,1,1,1,1,1,1,1,1,
				1,0,1,1,1,1,1,1,1,1,
				1,1,0,1,1,1,1,1,1,1,
				1,1,1,0,1,1,1,1,1,1,
				1,1,1,1,0,1,1,1,1,1,
				1,1,1,1,1,0,1,1,1,1,
				1,1,1,1,1,1,0,1,1,1,
				1,1,1,1,1,1,1,0,1,1,
				1,1,1,1,1,1,1,1,0,1,
				1,1,1,1,1,1,1,1,1,0,
				1,1,1,1,1,1,1,1,1,1)

	act_mat <- calcActivation(T_nw, b, n, K, flag_gen_data=TRUE)

	checkEquals(true_result, act_mat)
}
