#
# caculate activation matrix
#
calcActivation <- function(T_nw, b, n, K, flag_gen_data=FALSE) {

  kds <- matrix(b, nrow=n, ncol=K)
  activation_mat <- matrix(NA, nrow=n, ncol=K)
  
  for (k in 1:dim(kds)[2]) {
		nw <- T_nw
		inflow <- vector("list", length=n)
		in_deg <- apply(abs(nw), 2, sum)  # in-degree for each edge
		root_nodes <- which(in_deg == 0)  # find root_nodes
		
		# if no root nodes: stop
		if (length(root_nodes) == 0) {
			cat("Error: there are no root nodes\n")
		}
		done <- vector()
		
		# process root_nodes: then delete their outgoing edges: set new indegree
		for (i in root_nodes) {
			if (kds[i,k] != 0) {  # root nodes are active if not kd 
				inflow[[i]] <- c(inflow[[i]], 1)
				children <- which(nw[i, ] != 0)
				
				for (c in children) {
					if (kds[c,k] == 0) {  # if children kd -> they are inactive
						nw[c, ] <- 0
						done <- c(done, c)
					}
					else{
						inflow <- .setInflow_act(nw[i,c], c, inflow, flag_gen_data)
						in_deg[c] <- in_deg[c] - 1
					}
				}
			}
			nw[i, ] <- 0
			in_deg <- apply(abs(nw), 2, sum)
			done <- c(done, i)
		}
		
		# now proceed with nodes where in_deg is zero and which are not yet done
		ids_tmp <- which(in_deg == 0)
		ids <- ids_tmp[!ids_tmp %in% done]
		
		while (length(ids) > 0) {
			for (i in ids) {
				parents <- which(T_nw[, i] != 0)
				for (pa in parents) {
					if (sum(inflow[[pa]]) > 0) 
						inflow <- .setInflow_act(T_nw[pa,i], i, inflow, flag_gen_data)
				}
				children <- which(nw[i, ] != 0)
				for (c in children) {
					if (kds[c,k] == 0) {  # if children kd -> they are inactive
						nw[c,] <- 0
						in_deg <- apply(abs(nw), 2, sum)
						done <- c(done, c)
					}
					else{
						if (sum(inflow[[i]]) > 0)
							inflow <- .setInflow_act(nw[i,c], c, inflow, flag_gen_data)
						in_deg[c] <- in_deg[c]-1
					}
				}
				done <- c(done, i)
				ids_tmp <- which(in_deg == 0)
				ids <- ids_tmp[!ids_tmp %in% done]
			}
		}
		
		# if no indegree is zero and some nodes are still undone: there is a loop
		undone <- which(!seq(1, n) %in% done)
		incom <- unlist(lapply(inflow, sum))
		
		while (length(undone) > 0) {
			ids <- undone[incom[undone] != 0]  # if any undone node which is not root has inflow>0 start there
			if (length(ids) > 0) {
				for (i in ids) {
					if (kds[i,k] != 0) {  # if node is not kd: just let inflow like it is
						children <- which(nw[i, ] != 0)
						for (c in children) {
							if (kds[c,k] != 0) {  # if children kd -> they are inactive
								if (sum(inflow[[i]]) > 0) 
									inflow <- .setInflow_act(nw[i,c], c, inflow, flag_gen_data)
							}
						}
					}
					done <- c(done, i)
				}
			}
			else {  # if no node is active: rest is inactive
				for (i in undone) {
					inflow[[i]] <- c(inflow[[i]], 0)
					done <- c(done, i)
				}
			}
			undone <- which(!seq(1, n) %in% done)
			incom <- unlist(lapply(inflow, sum))
		}
		
		tmp <- unlist(lapply(inflow, sum))
		activation_mat[ ,k] <- apply(cbind(rep(0, n), tmp), 1, max)
	}
  activation_mat[activation_mat != 0] <- 1

  return(activation_mat)
}


.setInflow_act <- function(edge, node, inflow, flag_gen_data=FALSE) {
	
	if (edge > 0 | flag_gen_data) {
		inflow[[node]] <- c(inflow[[node]], 1)
	}
	else {
		inflow[[node]] <- c(inflow[[node]], -1)
	}
	
	return(inflow)
}

