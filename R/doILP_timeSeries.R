.doILP_timeSeries <- function(obs, delta, lambda, b, n, K, T_, annot, delta_type, prior=NULL, sourceNode=NULL,
															sinkNode=NULL, all.int=FALSE, all.pos=FALSE) {

	nConstr <- n*K*(T_-1)
	
  ## weight matrix of dim ((K*n)x(2n²+n)) (w_i^0)
  if(all.pos) {
		W <- matrix(0, nrow=nConstr, ncol=n*n+n)
	}
  else {
		W <- matrix(0, nrow=nConstr, ncol=2*n*n+n)
	}
  colnames(W) <- annot
  
  # direction of inequation
  f.dir <- rep("<=", nConstr)
  
  # Vector of numeric values for the right-hand sides of the constraints
  bvec <- rep(0, nConstr)
  J <- seq(1,n)
  count <- 1
	
	lp_problem_data <- list(W=W, bvec=bvec, f.dir=f.dir, count=count)
  if (delta_type == "perGene") {
		
		if (all.int)
			delta <- rep(1, n)
			
		for (t in 2:T_) {
			for (k in 1:K) {
				for (i in 1:n) {
					delta_i <- delta[i]
					lp_problem_data <- .loopThroughMatrix_ts(lp_problem_data, delta_type, i, k, t, n, 
																													obs, b, delta, delta_i, annot, J, all.pos)
				} 
			}
		} 
	}
	else if (delta_type == "perGeneExp") {
  
		if (all.int)
			delta <- matrix(rep(1, n*K), nrow=n, ncol=K)
		
		for (t in 2:T_) {
			for (k in 1:K) {
				for (i in 1:n) {
					delta_i <- delta[i,k]
					lp_problem_data <- .loopThroughMatrix_ts(lp_problem_data, delta_type, i, k, t, n, 
																													obs, b, delta, delta_i, annot, J, all.pos)
				} 
			} 
		} 
	}
	else if (delta_type == "perGeneTime") {
  
		if (all.int)
			delta <- matrix(rep(1, n*T_), nrow=n, ncol=T_)
		
		for (t in 2:T_) {
			for (k in 1:K) {
				for (i in 1:n) {
					delta_i <- delta[i,t]
					lp_problem_data <- .loopThroughMatrix_ts(lp_problem_data, delta_type, i, k, t, n, 
																													obs, b, delta, delta_i, annot, J, all.pos)
				} 
			} 
		}
	}
	else if (delta_type == "perGeneExpTime") {
  
		if (all.int)
			delta <- array(rep(1, n*K*T_), c(n,K,T_))
		
		for (t in 2:T_) {
			for (k in 1:K) {
				for (i in 1:n) {
					delta_i <- delta[i,k,t]
					lp_problem_data <- .loopThroughMatrix_ts(lp_problem_data, delta_type, i, k, t, n, 
																													obs, b, delta, delta_i, annot, J, all.pos)
				}
			} 
		} 
	}
  
  
	lp_problem_data <- .setSlackVariables_ts(lp_problem_data, n, nConstr, lambda, annot,all.pos)
	lp_problem_data <- .setPriors_ts(lp_problem_data, delta, n, annot, sourceNode, sinkNode, prior, all.pos)
	
  W <- lp_problem_data$W
  bvec <- lp_problem_data$bvec
  f.dir <- lp_problem_data$f.dir
  cvec <- lp_problem_data$cvec

  ## Maximize the gross margin
  res <- lp("min",cvec,W,f.dir,bvec,all.int=all.int) 
  ## min - direction of optimization
  ## cvec - objective function (Numeric vector of coefficients of objective function)
  ## W - Matrix of numeric constraint coefficients, one row per constraint, one column per variable
  ## f.dir vector of character strings giving the direction of the constraint
  ## bvec - vector of numeric values for the right-hand sides of the constraints

  return(res)
}


.loopThroughMatrix_ts <- function(lp_problem_data, delta_type, i, k, t, n, obs, 
																				 b, delta, delta_i, annot, J, all.pos) {

	bvec <- lp_problem_data$bvec
	f.dir <- lp_problem_data$f.dir
	count <- lp_problem_data$count
	
	if (b[(k-1)*n + i] == 1) {  # if the entry in b is 1 then the gene is active
		if (!is.na(obs[i,k,t])) {  # if the observation=NA, just do nothing
			
			if (obs[i,k,t]>= delta_i) {  # if observation of gene i after knockdown k is active
				if (all.pos) {  # set offset parameter (baseline of gene i)
					lp_problem_data$W[count, i+(n*n)] <- 1
					
					for (j in J[J != i]) {  # sum
						delta_j <- .setDeltaValue_ts(delta_type, delta, j, k, t)
						lp_problem_data <- .setPosMatrixEntries_ts(lp_problem_data, i, j, k, t, n, obs, b, delta_j, annot) 
					}
				}
				else {
					lp_problem_data$W[count,i+(2*n*n)] <- 1
					for (j in J[J != i]) {  # sum
						delta_j <- .setDeltaValue_ts(delta_type, delta, j, k, t)
						lp_problem_data <- .setPosNegMatrixEntries_ts(lp_problem_data, i, j, k, t, n, obs, b, delta_j, annot)
					} 
				}
				f.dir[count] <- ">="
				bvec[count] <- delta_i
			}
			
			if (obs[i,k,t] < delta_i) {  # if observation of gene i after knockdown k is NOT active
				if (all.pos) {
					lp_problem_data$W[count,i+(n*n)] <- 1  # set offset parameter (baseline of gene i)
					
					for (j in J[J != i]) {  # sum
						delta_j <- .setDeltaValue_ts(delta_type, delta, j, k, t)
						lp_problem_data <- .setPosMatrixEntries_ts(lp_problem_data, i, j, k, t, n, obs, b, delta_j, annot) 
					}
				}
				else {
					lp_problem_data$W[count,i+(2*n*n)] <- 1
					for (j in J[J != i]) {  # sum
							delta_j <- .setDeltaValue_ts(delta_type, delta, j, k, t)
							lp_problem_data <- .setPosNegMatrixEntries_ts(lp_problem_data, i, j, k, t, n, obs, b, delta_j, annot)
					}
				}
			f.dir[count] <- "<="
			bvec[count] <- 0
			}
		}
	}
	count <- count+1
	
	lp_problem_data$bvec <- bvec
	lp_problem_data$f.dir <- f.dir
	lp_problem_data$count <- count
	
	return(lp_problem_data)
}


.setDeltaValue_ts <- function(delta_type, delta, j, k, t) {
	
	if (delta_type == "perGene") {
		delta_j <- delta[j]
	}
	else if (delta_type == "perGeneExp") {
		delta_j <- delta[j,k]
	}
	else if (delta_type == "perGeneTime") {
		delta_j <- delta[j,t-1]
	}
	else if (delta_type == "perGeneExpTime") {
		delta_j <- delta[j,k,t-1]
	}
	
	return(delta_j)
}


.setPosMatrixEntries_ts <- function(lp_problem_data, i, j, k, t, n, obs, b, delta_j, annot) {

	W <- lp_problem_data$W
	count <- lp_problem_data$count
	
	idPos <- which(annot == paste("w+", j, i, sep="_"))
	if (!is.na(obs[j,k,t-1])) {
		if ((obs[j,k,t-1] >= delta_j) & (b[(k-1)*n + j] == 1)){
			W[count,idPos] <- obs[j,k,t-1]
		}
		else {
			W[count,idPos] <- 0
		}
	}
	else {
		W[count,idPos] <- NA
	}
	
	lp_problem_data$W <- W
	
	return(lp_problem_data)
}


.setPosNegMatrixEntries_ts <- function(lp_problem_data, i, j, k, t, n, obs, b, delta_j, annot) {

	W <- lp_problem_data$W
	count <- lp_problem_data$count
	
	#	 positive parameter
	idPos <- which(annot == paste("w+", j, i, sep="_"))
	idNeg <- which(annot == paste("w-", j, i, sep="_"))
	
	if (!is.na(obs[j,k,t-1])) {
		if ((obs[j,k,t-1] >= delta_j) & (b[(k-1)*n + j] == 1)){
			W[count,idPos] <- obs[j,k,t-1]
			W[count,idNeg] <- -obs[j,k,t-1]
		}
		else{
			W[count,idPos] <- 0
			W[count,idNeg] <- 0
		}
	}
	else{
		W[count,idPos] <- NA
		W[count,idNeg] <- NA
	}

	lp_problem_data$W <- W
	
	return(lp_problem_data)
}


.setSlackVariables_ts <- function(lp_problem_data, n, nConstr, lambda, annot, all.pos){

	W <- lp_problem_data$W
	f.dir <- lp_problem_data$f.dir

	if (lambda != 0) {
		sl <- matrix(0, nrow=nConstr, ncol=nConstr)
		annot_s <- paste("s",seq(1, nConstr), sep="_")
		colnames(sl) <- annot_s
		
		# attention: for the constr. where observation is smaller than threshold
		xi <- vector()
		for (j in 1:length(f.dir)) {
			if(f.dir[j] == ">=") 
				xi <- c(xi ,0)
			if(f.dir[j] == "<=") 
				xi <- c(xi, -1)
		}
		diag(sl) <- xi
		W <- cbind(W, sl)
		if (all.pos) {
			cvec <- c(rep(1, n*n), rep(1, n), rep(1/lambda, nConstr))
			# self-activation is not allowed
			id_self <- c(which(annot == paste("w+", seq(1, n), seq(1, n), sep="_")))
		}
		else {
			cvec <- c(rep(1, n*n), rep(1, n*n), rep(1, n), rep(1/lambda, nConstr))
			# self-activation is not allowed
			id_self <- c(which(annot == paste("w+", seq(1, n) ,seq(1, n), sep="_")),
									 which(annot == paste("w-", seq(1, n), seq(1, n), sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot, annot_s)
		
	}
	else {
		if (all.pos) {
			cvec <- c(rep(1, n*n), rep(1, n)) 
			# self-activation is not allowed
			id_self <- c(which(annot == paste("w+", seq(1, n), seq(1, n), sep="_")))
		}
		else {
			cvec <- c(rep(1, n*n), rep(1, n*n), rep(1, n)) 
			# self-activation is not allowed
			id_self <- c(which(annot == paste("w+", seq(1, n), seq(1, n), sep="_")),
									 which(annot == paste("w-", seq(1, n), seq(1, n), sep="_")))
		}
		cvec[id_self] <- 0
		names(cvec) <- c(annot)
  }

  lp_problem_data$W <- W
	lp_problem_data$f.dir <- f.dir
	lp_problem_data$cvec <- cvec
	
	return(lp_problem_data)
}


.setPriors_ts <- function(lp_problem_data, delta, n, annot, sourceNode, sinkNode, prior, all.pos){

	W <- lp_problem_data$W
	bvec <- lp_problem_data$bvec
	f.dir <- lp_problem_data$f.dir
	cvec <- lp_problem_data$cvec
	
  ## condition that each node which is not End hast at least delta[i] outgoing edges
  if (!is.null(sinkNode)) {
		W_tmp1 <- vector()
		gene_tmp <- seq(1, n)[-sinkNode]
		
		for (i in gene_tmp) {
			# outgoing edge can come from all nodes except itself
			tmp <- seq(1, n)[-i]
			
			if (length(tmp) > 1) {
				# for negative and positive parameter
				annot_pos <- paste("w+", i, tmp, sep="_")
				
				if (!all.pos) 
					annot_neg <- paste("w-", i, tmp, sep="_")
				
				add_row <- rep(0, length(cvec))
				add_row[which(annot %in% annot_pos)] <- 1
				
				if (!all.pos) 
					add_row[which(annot %in% annot_neg)] <- 1
				
				W_tmp1 <- rbind(W_tmp1, as.double(add_row))
				bvec <- c(bvec, delta[i])
				f.dir <- c(f.dir, ">=")
			}
		}
		W <- rbind(W, W_tmp1)
  }
  
  ## conditions that each node which is not Start has at least delta[i] incoming edges
  if (!is.null(sourceNode)) {
		W_tmp2 <- vector()
		gene_tmp <- seq(1, n)[-sourceNode]
		
		for (i in gene_tmp) {
			# incoming edge can come from all nodes except itself
			tmp <- seq(1, n)[-i]
			
			if (length(tmp) > 1) {
				annot_pos <- paste("w+", tmp, i, sep="_")
				if (!all.pos) 
					annot_neg <- paste("w-", tmp, i, sep="_")
				
				add_row <- rep(0, length(cvec))
				add_row[which(annot %in% annot_pos)] <- 1
				
				if (!all.pos) 
					add_row[which(annot %in% annot_neg)] <- 1
				
				W_tmp2 <- rbind(W_tmp2, as.double(add_row))
				bvec <- c(bvec, delta[i])
				f.dir <- c(f.dir, ">=")
			}
		}
		W <- rbind(W, W_tmp2)
  }

  ## if there is a prior
  if (!is.null(prior)) {
		for (i in 1:length(prior)) {
		
			tmp <- rep(0, dim(W)[2])
			tmp[which(prior[[i]][1] == annot)] <- as.double(prior[[i]][2])
			W <- rbind(W, tmp)
			bvec <- c(bvec, as.double(prior[[i]][4]))
			f.dir <- c(f.dir, prior[[i]][3])
		}
  }
  
  lp_problem_data$W <- W
	lp_problem_data$bvec <- bvec
	lp_problem_data$f.dir <- f.dir
	
	return(lp_problem_data)
}


