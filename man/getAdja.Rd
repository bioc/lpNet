\name{getAdja}
\alias{getAdja}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Adjacency Matrix.
}
\description{
The function returns the adjacency matrix of the network computed with the "doILP" function.
}
\usage{
getAdja(res, n, annot=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
Result returned by the "doILP" function.
}
  \item{n}{
Integer: the number of nodes of the inferred network.
}
  \item{annot}{
Vector of character strings: the annotation of the edges as returned by "getEdgeAnnot".
}
}

\value{
Numeric matrix: the adjacency matrix of the network.
}

\seealso{
\code{\link{doILP}}
}
\examples{
n <- 3 # number of genes
K <- 4 # number of experiments
T_ <- 4 # number of time points

# generate random observation matrix
obs <- array(rnorm(n*K*T_), c(n,K,T_))

baseline <- c(0.75, 0, 0)

delta <- rep(0.75, n)

# perturbation vector, entry is 0 if gene is inactivated and 1 otherwise
b <- c(0,1,1, # perturbation exp1: gene 1 perturbed, gene 2,3 unperturbed
       1,0,1, # perturbation exp2: gene 2 perturbed, gene 1,3 unperturbed
       1,1,0, # perturbation exp3....
       1,1,1)
            
delta_type <- "perGene"
lambda <- 1/10
annot <- getEdgeAnnot(n)

#infer the network
res <- doILP(obs, delta, lambda, b, n, K, T_, annot, delta_type, prior=NULL,
             sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE, flag_time_series=TRUE)

# make the adjacency matrix
adja <- getAdja(res, n)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{adjacency }


