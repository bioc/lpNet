\name{generateTimeSeriesNetStates}
\alias{generateTimeSeriesNetStates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Time Series Network States
}
\description{
The function returns all gene states for each network state in time-series data. The signalling propagates downstream one edge per time-point. The stopping criteria is when all edges have been active at least once, so that infinite loops are avoided. The number of time points for the data can be defined by the user or not, if not the number of time points will be the same as the number of different network states. If the number of time points is defined by the user, network states will be either repeated or removed, so that there are as many network states as time points.
}
\usage{
generateTimeSeriesNetStates(nw_und, b, n, K, T_user=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nw_und}{
Numeric matrix: the adjacency matrix representing the underlying network.
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
  \item{T_user}{
Integer definining the number of time points in the network.
}
}

\value{
List containing an array with all nodes states and the number of time points.
}

\examples{
n <- 3 # number of genes
K <- 4 # number of experiments

# perturbation vector, entry is 0 if gene is inactivated and 1 otherwise
b <- c(0,1,1, # perturbation exp1: gene 1 perturbed, gene 2,3 unperturbed
       1,0,1, # perturbation exp2: gene 2 perturbed, gene 1,3 unperturbed
       1,1,0, # perturbation exp3....
       1,1,1)

# adjacency matrix       
nw_und <- matrix(c(0,1,0,
                   0,0,1,
                   0,0,0), nrow=n, ncol=n, byrow=TRUE)

generateTimeSeriesNetStates(nw_und,b, n, K, T_user=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{matrix summarization }


