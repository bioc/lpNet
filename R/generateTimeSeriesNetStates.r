#
# function that actually builds the time series data (node states at 
# each time point) with the desired number if time points
#
generateTimeSeriesNetStates <- function(nw_und, b, n, K, T_user=NULL) {
	
	# get the number of steps required for the signal to propagate 
	# from the source to the end nodes
	
	T_ <- .getNTimePoints(nw_und, n, 1)
	
	# get the activation matrix that doesn't take into account the 
	# incoming edges sign
	act_mat <- calcActivation(nw_und, b, n, K, flag_gen_data=TRUE)
	
	# at t <- {0,1}, no edges are active
	active_nw <- list()
	active_nw_temp <- matrix(0, nrow=n, ncol=n)
	active_nw[[1]] <- active_nw_temp
	active_nw[[2]] <- active_nw_temp
	
	# at t <- 0 no nodes are active
	node_stat <- list()
	active_nodes <- rep(0, n)
	node_stat[[1]] <- matrix(rep(active_nodes, K), nrow=n, ncol=K)
	
	# initialize matrix for node states
	for (t in 2:T_) 
		node_stat[[t]] <- matrix(NA, nrow=n, ncol=K)
	
	# generate the node states for each experiment
	for (k in 1:K) {
		# get the parent node, 
		in_deg <- apply(abs(nw_und), 2, sum)
		parent_nodes <- which(in_deg == 0)
		silenced_parents <- which(parent_nodes %in% which(act_mat[ ,k]==0))
		
		# the silenced parent nodes cannot influence any other nodes
		if (length(silenced_parents) >= 1) {
			parent_nodes <- parent_nodes[-silenced_parents] 
	 }
		root_nodes <- parent_nodes
		non_root_nodes <- which(!(seq(1, n) %in% root_nodes))
		
		active_nw_temp <- matrix(0, nrow=n, ncol=n)
		
		# at t <- 1 only non silenced root nodes are active
		active_nodes <- rep(0, n)
		active_nodes[root_nodes] <- 1
		node_stat[[2]][ ,k] <- active_nodes
		
		# initialize variables
		count <- 1
		visited_nodes <- vector()
		edges_added <- vector()
		n_edges <- length(which(nw_und!=0))
		
		for (t in 3:T_) {
		
			all_children <- c()

			for (i in parent_nodes) {
				children <- c()
				poss_child <- which(nw_und[i, ] != 0)

				# condition to avoid loops in net construction
				# if the edge i -> j is not zero, then j is not  child anymore
				for (j in poss_child) {
					if (active_nw[[1+count]][i,j]==0){
						children <- c(children, j)
					}
				}
				all_children <- c(all_children, children)
				active_nw_temp[i, children] <- nw_und[i,children]
			}

			active_nw[[2+count]] <- active_nw_temp
			active_nodes[all_children] <- 1
		
			# set nodes as active or inactive according to their incoming edges
			# if the sum of incoming edges is positive, node is active
			# if the sum of incoming edges is negative, node is inactive
			# the active/inactive state of a root node doesn't 
			# change with net evolution, as it has no incoming edges
			for (i in non_root_nodes) {
				if ((sum(active_nw_temp[ ,i]) > 0) & (act_mat[i,k] == 1)) {
					active_nodes[i] <- 1
				}
				else {
					active_nodes[i] <- 0
				}
			}
			# if a given node is inactive, then it has no active outgoing edges
			for (i in non_root_nodes) {
				if (active_nodes[i] == 0) 
					active_nw_temp[i, ] <- 0
			}

			node_stat[[2+count]][ ,k] <- active_nodes

			parent_nodes <- which(active_nodes > 0)
			count <- count + 1
			

		} # end of T_
	} # end of k	
	
	if (!is.null(T_user)) {
		# build final array with node states
		node_state_vec <- array(NA, c(n,K,T_user))
	
		# if the number of time points is greater than the total number of 
		# steps it takes for the signal to propagate to the sink nodes: 
		# randomly repeat some node state vectors
		if (T_user > T_) {
			to_repeat <- sample(seq(1, T_), T_user-T_, replace=TRUE)
			time_points <- sort(c(seq(1, T_), to_repeat))
			i <- 1
			for (t in time_points) {
				node_state_vec[ , ,i] <- node_stat[[t]]
				i <- i + 1
			}
		}
		
		# if the number of time points is less than the total number of 
		# steps it takes for the signal to propagate to the sink nodes: 
		# randomly delete some node state vectors
		else if (T_user < T_) {
			to_remove <- sample(seq(1, T_), T_-T_user, replace=TRUE)
			time_points <- seq(1 ,T_)
			time_points <- time_points[-to_remove]
			i <- 1
			for (t in time_points) {
				node_state_vec[ , ,i] <- node_stat[[t]]
				i <- i + 1
			}
		}
	}
	else {
		# build final array with node states
		node_state_vec <- array(NA, c(n,K,T_))
		for (t in 1:T_)
			node_state_vec[ , ,t] <- node_stat[[t]]
	}
	
	return(list(node_state_vec=node_state_vec, T_=dim(node_state_vec)[3]))
}


#
# this function gets the total number of steps it takes for the signal to
# propagate from the root to the root to the sink nodes, when no nodes
# are silenced
#
.getNTimePoints <- function(nw_und, n, K) {
	
	# active_nw contains only the active connections at each time point
	# at t <- {0,1} there are no active connections
	active_nw <- list()
	active_nw_temp <- matrix(0, nrow=n, ncol=n)
	active_nw[[1]] <- active_nw_temp
	active_nw[[2]] <- active_nw_temp
	
	# at t <- 0 no nodes are active
	node_stat <- list()
	active_nodes <- rep(0, n)
	node_stat[[1]] <- matrix(rep(active_nodes, K), nrow=n, ncol=K)
	
	# get the root, non root, and the temporary parent nodes
	in_deg <- apply(abs(nw_und), 2, sum)
	parent_nodes <- which(in_deg == 0)
	root_nodes <- parent_nodes
	non_root_nodes <- which(!(seq(1, n) %in% root_nodes))
	
	# at t <- 1, the root nodes are active
	active_nodes <- rep(0, n)
	active_nodes[root_nodes] <- 1
	node_stat[[2]] <- matrix(rep(active_nodes, K), nrow=n, ncol=K)
	
	# initialize variables
	count <- 1
	visited_nodes <- vector()
	edges_added <- vector()
	n_edges <- length(which(nw_und != 0))
	
	# the loop ends when all edges have been added or when there are no more children
	while ((length(unique(edges_added)) < n_edges)) {
		all_children <- c()

		for (i in parent_nodes) {
			children <- c()
			poss_child <- which(nw_und[i, ] != 0)

			# condition to avoid loops in net construction:
			# if the edge i -> j is not zero, then j is not a child anymore
			for (j in poss_child) {
				if (active_nw[[count+1]][i,j] == 0) {
					children <- c(children, j)
				}
			}
			all_children <- c(all_children, children)
			active_nw_temp[i,children] <- nw_und[i,children]
			edges_added <- c(edges_added, sprintf("%s->%s", i, children))
		}
		
		# if there are no more children break
		if (length(all_children) == 0) {
			break
		}

		active_nw[[2+count]] <- active_nw_temp
		active_nodes[all_children] <- 1
	
		# set nodes as active or inactive according to their incoming edges
		# if the sum of incoming edges is positive, node is active
		# if the sum of incoming edges is negative, node is inactive
		for (i in non_root_nodes) {
			if ( sum(active_nw_temp[ ,i]) > 0) {
				active_nodes[i] <- 1
			}
			else {
				active_nodes[i] <- 0
			}
		}
		# if a given node is inactive, then it has no active outgoing edges
		for (i in non_root_nodes) {
			if (active_nodes[i] == 0)
				active_nw_temp[i, ] <- 0
		}

		# create the matrix with node states at t
		node_stat[[2+count]] <- matrix(rep(active_nodes, K), nrow=n, ncol=K)

		# the parent nodes are the active nodes at t
		parent_nodes <- which(active_nodes > 0)
		count <- count + 1
		
		edges_added <- unique(edges_added)
	}
	
	return(length(node_stat))
}

