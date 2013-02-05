\name{doILP}
\alias{doILP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Do The Network Inference With The Linear Programming Approach.
}
\description{
This function converts observation data into a linear programming problem.
}
\usage{
doILP(obs, delta, lambda, b, n, K, annot, 
  prior = NULL, sourceNode=NULL, sinkNode=NULL, all.int = FALSE, all.pos = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
Numeric matrix: the given observations.
}
  \item{delta}{
Numeric vector defining the thresholds for each gene to determine its observation to be active or inactive.
}
  \item{lambda}{
Numeric value defining the penalty parameter lambda. It can range from zero to infinity and it controls the introduction of slack variables in the network inference lp model.
}
  \item{b}{
Vector of 0/1 values describing the experiments (entry is 0 if gene is inactivated in the respetive experiment and 1 otherwise). The measurements of the genes of each experiment are appended as a long vector.
}
  \item{n}{
Integer: number of genes.
}
  \item{K}{
Integer: number of perturbation experiments.
}
  \item{annot}{
Vector of character strings: the annotation of the edges as returned by "getEdgeAnnot".
}
  \item{prior}{
Prior knowledge, given as a list of constraints. Each constraint consists of a vector with four entries describing the prior knowledge of one edge. For example the edge between node 1 and 2, called w+_1_2, is defined to be bigger than 1 with constraint c("w+_1_2",1,">",2). The first entry specifies the annotation of the edge (see function "getEdgeAnnot") and the second defines the coefficient of the objective function (see parameter "objective.in" in the "lp" function of the package "lpSolve"). Furthermore, the third, respectively the fourth elements give the direction, respectively the right-hand side of the constraint (see the parameters "const.dir", respectively "const.rhs" in the "lp" function of the package "lpSolve").
}
  \item{sourceNode}{
Integer vector: indices of the known source nodes.
}
  \item{sinkNode}{
Integer vector: indices of the known sink nodes.
}
  \item{all.int}{
Logical: should all variables be integer? Corresponds to an Integer Linear Program (see "lp" function in package "lpSolve"). Default: FALSE.
}
  \item{all.pos}{
Logical: should all variables be positive? Corresponds to learning only activating edges. Default: FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An lp object. See "lp.object" in package "lpSolve" for details.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%\code{\link{lp}}
}
\examples{
n <- 5 # number of genes
K <- 7 # number of perturbation experiments
annot <- getEdgeAnnot(n)

# generation of random observation matrix
obs <- matrix(rnorm(35),nrow=n,ncol=K)
# define delta 
delta <- apply(obs,1,mean,na.rm=TRUE)

# activation (knockdown) vector (entry is 0 if gene is inactivated in the respetive experiment)
b <- c(0,1,1,1,1,
1,0,1,1,1,
1,1,0,1,1,
1,1,1,0,1,
1,1,1,1,0,
1,0,0,1,1,
1,1,1,1,1)

res <- doILP(obs,delta,lambda=1,b,n,K,annot)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ network inference }
\keyword{ linear programming approach }% __ONLY ONE__ keyword per line
