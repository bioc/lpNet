test.getSampleAdjaMAD <- function() {
	
	n <- 3
	K <- 4
	annot <- getEdgeAnnot(n)
	annot_node = seq(1,n)
	
	true_result <- matrix(c(0, 0.7947368, 0.0000000, 
													0, 0.0000000, 0.0000000, 
													0, 0.0000000, 0.0000000), nrow=n, ncol=n, byrow=TRUE)
	colnames(true_result) <- rownames(true_result) <- annot_node
	
	edges_all <- matrix(c(0.7947368, 0.7947368, 0, 0.0000000, 0, 0.0000000,
												0.0000000, -1.1411606, 0, 1.9358974, 0, 0.0000000,
												0.0000000, -1.1411606, 0, 1.9358974, 0, 1.3482143,
												0.7947368, 0.7947368, 0, 0.0000000, 0, 0.0000000,
												0.7947368, 0.0000000, 0, 0.7947368, 0, 0.0000000,
												0.7947368, 0.7947368, 0, 0.0000000, 0, 0.0000000,
												-0.5534774, -1.1411606, 0, 1.9358974, 0, 1.3482143,
												0.7947368, -1.1411606, 0, 1.9358974, 0, 0.0000000,
												0.7947368, -1.1411606, 0, 1.9358974, 0, 0.0000000,
												0.3262604, -0.7947368, 0, 0.7947368, 0, 0.7947368,
												1.9358974, 0.0000000, 0, -1.3482143, 0, -1.9358974,
												1.9358974, 0.0000000, 0, 0.0000000, 0, -1.9358974), nrow=n*K, ncol=n*(n-1), byrow=TRUE)
	
	colnames(edges_all) <- c("1->2", "1->3", "2->1", "2->3", "3->1", "3->2")
	
	sampleAdjaMAD = getSampleAdjaMAD(edges_all, n, annot_node, method=median, method2=mad, septype="->")

	checkEquals(true_result, sampleAdjaMAD, tolerance=0.00001)
}
