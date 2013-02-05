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
getAdja(result, numnodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{result}{
Result returned by the "doILP" function.
}
  \item{numnodes}{
Integer: the number of nodes of the inferred network.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric matrix: the adjacency matrix of the network.
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
\code{\link{doILP}}
}
\examples{
# generate random observation matrix with 5 experiments and 5 genes
obs <- matrix(rnorm(5*5,1,0.1),nrow=5,ncol=5)
# define the perturbations
b <- c(0,1,1,1,1,
1,0,1,1,1,
1,1,0,1,1,
1,1,1,0,1)
n <- 5 # number of genes
K <- 4 # number of knockdowns
# annotation of the edges
annot <- getEdgeAnnot(n)
# define delta
delta <- rep(1,5)
# infer the network
res <- doILP(obs,delta,lambda=1,b,n,K,annot)

# make the adjacency matrix
adja <- getAdja(res,n)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{adjacency }


