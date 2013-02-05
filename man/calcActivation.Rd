\name{calcActivation}
\alias{calcActivation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Activation Matrix
}
\description{
Calculate the activation matrix assuming that the signaling is deterministically propagated along the network. For a given network and given perturbation experiments the theoretical states of the genes are computed. So if a gene has been silenced in an experiment, then the state of this gene is assumed to be inactive, otherwise it is active, if its inflow (coming from parent nodes) is activating. Cycles within a network are not resolved, therefore this function can be used only for networks without cycles.
}
\usage{
calcActivation(T_nw, b, n, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T_nw}{
Adjacency matrix: the network which is used to compute the activities and inactivities.
}
  \item{b}{
Vector of 0/1 values describing the experiments (entry is 0 if gene is inactivated in the respetive experiment and 1 otherwise). The measurements of the genes of each experiment are appended as a long vector.
}
  \item{n}{
Integer: number of genes.
}
  \item{K}{
Integer: number of pertrubations experiments.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Matrix of 0/1 values; rows corresponding to genes, columns to experiments. If an entry is 1, it means that the corresponding gene is active in the corresponing experiment and inactive otherwise.
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
n <- 5 # number of genes
K <- 7 # number of perturbations experiments

# perturbation vector, entry is 0 if gene is inactivated and 1 otherwise
b <- c(0,1,1,1,1, # perturbation exp1: gene 1 perturbed, gene 2-5 unperturbed
1,0,1,1,1,        # perturbation exp2: gene 2 perturbed, gene 1,3,4,5 unperturbed
1,1,0,1,1,        # perturbation exp3....
1,1,1,0,1,
1,1,1,1,0,
1,0,0,1,1,
1,1,1,1,1)
# example network
T_nw <- matrix(c(0,1,1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0),
  nrow=n,ncol=n,byrow=TRUE)
# compute theoretical activation of genes from example network with given perturbations
act_mat <- calcActivation(T_nw,b,n,K)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{activation }

