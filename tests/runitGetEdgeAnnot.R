test.getEdgeAnnot <- function() {
	
	true_result = c("w+_1_1", "w+_1_2", "w+_1_3", "w+_2_1", "w+_2_2", "w+_2_3", "w+_3_1", "w+_3_2", "w+_3_3",
									"w-_1_1", "w-_1_2", "w-_1_3", "w-_2_1", "w-_2_2", "w-_2_3", "w-_3_1", "w-_3_2", "w-_3_3",
									"w_1_^_0", "w_2_^_0", "w_3_^_0")
	
	n <- 3
	edge_annot <- getEdgeAnnot(n, allpos=FALSE)

	checkEquals(true_result, edge_annot)
}


test.getEdgeAnnotAllPos <- function() {
	
	true_result = c("w+_1_1", "w+_1_2", "w+_1_3", "w+_2_1", "w+_2_2", "w+_2_3", "w+_3_1", "w+_3_2", "w+_3_3",
									"w_1_^_0", "w_2_^_0", "w_3_^_0")
	
	n <- 3
	edge_annot <- getEdgeAnnot(n, allpos=TRUE)

	checkEquals(true_result, edge_annot)
}

