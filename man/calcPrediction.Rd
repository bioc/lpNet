\name{calcPrediction}
\alias{calcPredictionLOOCV}
\alias{calcPredictionKfoldCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Predicted Observation.
}
\description{
Calculate the predicted observation of a perturbation experiment. If observations of an experiment are missing this function can be used to determine for a given network the predicted outcome. The missing measurement is predicted from two normal distributions, one for observations coming from active and one coming from inactive genes. 
}
\usage{
calcPredictionLOOCV(kds, adja, obs, delta, rem_kd, rem_gene, 
  active_mu, active_sd, inactive_mu, inactive_sd)
calcPredictionKfoldCV(adja, b, n, K, 
  active_mu, active_sd, inactive_mu, inactive_sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kds}{
  Matrix of 0/1 values representing the perturbation experiments (entry is 0 if gene is inactivated in the respective experiment and 1 otherwise). One column corresponds to one experiment and each row represents a gene.
}
  \item{adja}{
Numeric matrix: the adjacency matrix of the given network.
}
  \item{obs}{
Numeric matrix: the observation matrix.
}
  \item{delta}{
Numeric vector defining the thresholds for each gene to determine its observation to be active or inactive.
}
  \item{rem_kd}{
Integer: the index of the experiment which is missing.
}
  \item{rem_gene}{
Integer: the index of the gene which is missing. Together with rem_kd, the two indices correspond to the observation of the gene rem_gene at experiment rem_kd.
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
  \item{b}{
Vector of 0/1 values describing the experiments (entry is 0 if gene is inactivated in the respetive experiment and 1 otherwise). The measurements of the genes of each experiment are appended as a long vector.
}
  \item{n}{
Integer: number of genes.
}
  \item{K}{
Integer: number of perturbations experiments.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric: the predicted observation(s).
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
\code{\link{kfoldCV}}, \code{\link{loocv}}
}
\examples{
n <- 5 # number of genes
K <- 7 # number of experiments

# perturbation vector, entry is 0 if gene is inactivated and 1 otherwise
b <- c(0,1,1,1,1, # perturbation exp1: gene 1 perturbed, gene 2-5 unperturbed
1,0,1,1,1,        # perturbation exp2: gene 2 perturbed, gene 1,3,4,5 unperturbed
1,1,0,1,1,        # perturbation exp3....
1,1,1,0,1,
1,1,1,1,0,
1,0,0,1,1,
1,1,1,1,1)
T_nw <- matrix(c(0,1,1,0,0,0,0,NA,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0),
  nrow=n,ncol=n,byrow=TRUE)

# define the parameters for the observation generated from the normal distributions
active_mu <- 0.9
inactive_mu <- 0.5
active_sd <- inactive_sd <- 0.01

# compute the predicted observation matrix for the "kfoldCV"
calcPredictionKfoldCV(T_nw, b, n, K, active_mu, active_sd, inactive_mu, inactive_sd)

#### LOOCV
# perturbation experiments written as a matrix
kds <- matrix(b,nrow=n,ncol=K)

# generate random observation matrix
obs <- matrix(rnorm(35),nrow=n,ncol=K)
# define delta
delta <- rep(0.75,5)
# define the observation to be removed
rem_kd <- 3
rem_gene <- 2
obs[rem_gene,rem_kd] <- NA
# compute the predicted value
calcPredictionLOOCV(kds, T_nw, obs, delta, rem_kd, rem_gene, 
  active_mu, active_sd, inactive_mu, inactive_sd)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mean squared error }

