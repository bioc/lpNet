\name{calcRangeLambda}
\alias{calcRangeLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Range Of Penalty Parameter Lambda.
}
\description{
The penalty parameter lambda can range from zero to infinity and it controls the introduction of slack variables in the network inference lp model. To limit the introduction of slack variables we restrict lambda to be not larger than lambdaMax (=the number of slack variables times the variance of all measurements given). This function computes the range from zero to lambdaMax with a given stepsize that increases exponentially.
}
\usage{
calcRangeLambda(obs, delta, delta_type, flag_time_series=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
Numeric matrix/array: the observation matrix/array. It can have up to 3 dimensions, where dimension 1 are the network nodes, dimension 2 are the perturbation experiments, and dimension 3 are the time points (if considered).
}
  \item{delta}{
Numeric: defines the thresholds for each gene to determine its observation to be active or inactive. This parameter can be either a numeric vector, a matrix, or a 3D array, depending on the specified delta_type.
}
  \item{delta_type}{
Character: can have the following values and meanings:
\itemize{
    \item "perGene" - the value of delta depends on the gene;
    \item "perGeneExp" - the value of delta depends on the gene and perturbation experiment;
    \item "perGeneTime" - the value of delta depends on the gene and time point;
    \item "perGeneExpTime" - the value of delta depends on the gene, perturbation experiment, and time point;
  }
}
  \item{flag_time_series}{
  Logical: specifies whether steady-state (FALSE) or time series data (TRUE) is used.  
}
}

\value{
Numeric vector of possible values for lambda.
}


\examples{
# generate random observation matrix with 5 experiments and 5 genes
obs <- matrix(rnorm(5*5, 1, 0.1), nrow=5, ncol=5)

# define delta to be 1 for each gene
delta <- rep(1, 5)
delta_type <- "perGene"

lambda_values <- calcRangeLambda(obs, delta, delta_type)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ penalty parameter }

