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
doILP(obs, delta, lambda, b, n, K, T_=NULL, annot, delta_type, prior=NULL,sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE, flag_time_series=FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
Numeric matrix/array: the given observation matrix/array. It can have up to 3 dimensions, where dimension 1 are the network nodes, dimension 2 are the perturbation experiments, and dimension 3 are the time points (if considered).
}
  \item{delta}{
Numeric: defining the thresholds for each gene to determine its observation to be active or inactive. This parameter can be either a numeric vector, a matrix, or a 3D array, depending on the specified delta_type.
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
  \item{T_}{
Integer: number of time points.
}
  \item{annot}{
Vector of character strings: the annotation of the edges as returned by "getEdgeAnnot".
}
  \item{delta_type}{
Character: can have the following values and meanings:
  - "perGene" - the value of delta depends on the gene;
  - "perGeneExp" - the value of delta depends on the gene and perturbation experiment;
  - "perGeneTime" - the value of delta depends on the gene and time point;
  - "perGeneExpTime" - the value of delta depends on the gene, perturbation experiment, and time point;
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
  \item{flag_time_series}{
  Logical: specifies whether steady-state (FALSE) or time series data (TRUE) is used.  
}
}

\value{
An lp object. See "lp.object" in package "lpSolve" for details.
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

res <- doILP(obs, delta, lambda, b, n, K, T_, annot, delta_type, prior=NULL,
             sourceNode=NULL, sinkNode=NULL, all.int=FALSE, all.pos=FALSE, flag_time_series=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ network inference }
\keyword{ linear programming approach }% __ONLY ONE__ keyword per line
