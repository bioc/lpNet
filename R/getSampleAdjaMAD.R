#
# sum up all inferred networks into a single one using method.
#  only edges whose "method2" (eg. mad) is greater than its "method" 
#  (e.g. median) will be included in the final network
#
getSampleAdjaMAD <- function(edges_all, n, annot_node, method=median, method2=mad, septype="->") {

  edge_med <- apply(edges_all, 2, method, na.rm=T)
  edge_mad <- apply(edges_all, 2, method2, na.rm=T)
  
  sample <- matrix(0, nrow=n, ncol=n)
  colnames(sample) <- rownames(sample) <- annot_node
  edgelist <- strsplit(colnames(edges_all), septype, fixed=T)
  
  for (i in 1:length(edge_med)) {
		id1 <- which(edgelist[[i]][1] == rownames(sample))
		id2 <- which(edgelist[[i]][2] == colnames(sample))
			
			if (abs(edge_med[i]) > abs(edge_mad[i]))
				sample[id1,id2] <- edge_med[i]
  }
  
  return(sample)
}



