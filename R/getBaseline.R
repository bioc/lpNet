#
# get baseline values when solving LP problem with lpSolve
#
getBaseline = function(res, n, allpos=FALSE) {

	if (allpos == FALSE) 
		baseline = res$solution[(2*n*n+1):(2*n*n+n)]
	else 
		baseline = res$solution[(n*n+1):(n*n+n)]
	
	return(baseline)
}
