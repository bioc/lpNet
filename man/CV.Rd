\name{CV}
\alias{kfoldCV}
\alias{loocv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation
}
\description{
Performs a stratified k-fold cross-validation or a Leave-One-Out cross-validation.
}
\usage{
kfoldCV(times, obs, n, b, K, delta, lambda, annot, annot_node, kfold,
  active_mu, active_sd, inactive_mu, inactive_sd,
  prior = NULL, sourceNode = NULL, sinkNode = NULL, allint = FALSE, allpos = FALSE)
loocv(times, obs, n, b, K, delta, lambda, annot, annot_node, 
  active_mu, active_sd, inactive_mu, inactive_sd, 
  prior = NULL, sourceNode = NULL, sinkNode = NULL, allint = FALSE, allpos = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
Integer: the number of times the cross-validation shall be performed.
}
  \item{obs}{
Numeric matrix: the measured observations.
}
  \item{n}{
Integer: number of genes.
}
  \item{b}{
Vector of 0/1 values describing the experiments (entry is 0 if gene is inactivated in the respetive experiment and 1 otherwise). The measurements of the genes of each experiment are appended as a long vector.
}
  \item{K}{
Integer: number of perturbation experiments.
}
  \item{delta}{
Numeric vector defining the thresholds to determine an observation active or inactive.
}
  \item{lambda}{
Numeric value defining the penalty parameter lambda. It can range from zero to infinity and it controls the introduction of slack variables in the network inference lp model. 
}
  \item{annot}{
Vector of character strings: the annotation of the edges as returned by "getEdgeAnnot".
}
  \item{annot_node}{
Vector of character strings: the annoation of the nodes.
}
  \item{kfold}{
Integer value of the number "k" in the k-fold cross-calidation.
}
  \item{active_mu}{
Numeric: the average value assumed for observations coming from activated nodes. The parameter active_mu and active_sd are used for predicting the observations of the normal distribution of activate genes.
}
  \item{active_sd}{
Numeric: the variation assumed for observations coming from activated nodes. The parameter active_mu and active_sd are used for predicting the observations of the normal distribution of activate genes.
}
  \item{inactive_mu}{
Numeric: the average value assumed for observations coming from inactivated nodes. The parameter inactive_mu and inactive_sd are used for predicting the observations of the normal distribution of inactivate genes.
}
  \item{inactive_sd}{
  Numeric: the variation assumed for observations coming from inactivated nodes. The parameter inactive_mu and inactive_sd are used for predicting the observations of the normal distribution of inactivate genes.
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
  \item{allint}{
Logical: should all variables be integer? Corresponds to an Integer Linear Program (see "lp" function in package "lpSolve"). Default: FALSE. 
}
  \item{allpos}{
Logical: should all variables be positive? Corresponds to learning only activating edges. Default: FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of 
 \item{MSE }{The mean squared error (MSE) of predicted and observed measurements of the corresponding cross-validation step.}
 \item{edges_all }{The learned edge weights for each cross-validation step.}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
annot_node <- seq(1,5)
n <- 5 # number of genes
K <- 7 # number of perturbation experiments

annot <- getEdgeAnnot(n)

# activation (knockdown) vector (entry is 0 if gene is inactivated in the respetive experiment)
b <- c(0,1,1,1,1,
1,0,1,1,1,
1,1,0,1,1,
1,1,1,0,1,
1,1,1,1,0,
1,0,0,1,1,
1,1,1,1,1)

T_nw <- matrix(c(0,1,1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0),nrow=n,ncol=n,byrow=TRUE)
colnames(T_nw) <- rownames(T_nw) <- annot_node

## calculate observation matrix with given parameters for the 
# Gaussian distributions for activation and deactivation
active_mu <- 0.95
inactive_mu <- 0.56
active_sd <- 0.1
inactive_sd <- 0.1
# use three replicates
replnum <- 3

# generation of random observation matrix
obs <- matrix(rnorm(35),nrow=n,ncol=K)

times <- kfold <- 10   # can be increased i.e. to 1000 to produce stable results

# define delta 
delta <- apply(obs,1,mean,na.rm=TRUE)

#### LOOCV
loocv_res <- loocv(times,obs,n,b,K,delta,lambda=1,annot,annot_node,
  active_mu,active_sd,inactive_mu,inactive_sd)

#### K-fold CV
kcv_res <- kfoldCV(times,obs,n,b,K,delta,lambda=1,annot,annot_node,kfold,
  active_mu,active_sd,inactive_mu,inactive_sd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cross-validation }
