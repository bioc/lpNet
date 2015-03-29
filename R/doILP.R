#
# build and solve LP problem
#
doILP <- function(obs, delta, lambda, b, n, K, T_=NULL, annot, delta_type, prior=NULL, 
									sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE, 
									flag_ss_v2=FALSE, flag_time_series=FALSE) {

	if (flag_time_series == FALSE) {
		if (flag_ss_v2 == TRUE) {
			res <- .doILP_steadyStateV2(obs, delta, lambda, b, n, K, T_=NULL, annot, delta_type, prior, 
																	sourceNode, sinkNode, all.int, all.pos)
		}
		else {
			res <- .doILP_steadyState(obs, delta, lambda, b, n, K, T_=NULL, annot, delta_type, prior, 
																sourceNode, sinkNode, all.int, all.pos)
		}
	}
	else {
		res <- .doILP_timeSeries(obs, delta, lambda, b, n, K, T_, annot, delta_type, prior, sourceNode,
														 sinkNode, all.int, all.pos)
	}
	
	return(res)
}
