\name{calcRangeLambda}
\alias{calcRangeLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Range Of Penalty Parameter Lambda.
}
\description{
The penalty parameter lambda can range from zero to infinity and it controls the introduction of slack variables in the network inference lp model. To limit the introduction of slack variables we restrict lambda to be not larger than lambdaMax (=the number of slack variables times the variance of all measurements given). This function computes the range from zero to lambdaMax with a given stepsize.
}
\usage{
calcRangeLambda(delta, obs, stepsize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
Numeric vector defining the thresholds for each gene to determine its observation to be active or inactive.
}
  \item{obs}{
Numeric matrix of the observed measurements.
}
  \item{stepsize}{
Numeric: giving the step size.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric vector of possible values for lambda.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%B. Knapp
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# generate random observation matrix with 5 experiments and 5 genes
obs <- matrix(rnorm(5*5,1,0.1),nrow=5,ncol=5)
# define delta to be 1 for each gene
delta <- rep(1,5)
lambda <- calcRangeLambda(delta,obs,stepsize=0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ penalty parameter }

