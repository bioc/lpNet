\name{getObsMat}
\alias{getObsMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Observation Matrix.
}
\description{
The function generates from the activation matrix, computed with "calcActivation", the observation matrix where active/inactive observations are generated from a normal distribution with the average and variation as given in the parameters.
}
\usage{
getObsMat(act_mat, active_mu, active_sd, inactive_mu, inactive_sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{act_mat}{
Matrix of 0/1 values called the activation matrix. Rows correspond to genes, columns to experiments. If an entry is 1, it means that the corresponding gene is active in the corresponing experiment and inactive otherwise.
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Numeric matrix: the observation matrix
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
\code{\link{calcActivation}}
}
\examples{
n <- 5 # number of genes
K <- 7 # number of knockdowns

# perturbation vector, entry is 0 if gene is inactivated and 1 otherwise
b <- c(0,1,1,1,1, # perturbation exp1: gene 1 perturbed, gene 2-5 unperturbed
1,0,1,1,1,        # perturbation exp2: gene 2 perturbed, gene 1,3,4,5 unperturbed
1,1,0,1,1,        # perturbation exp3....
1,1,1,0,1,
1,1,1,1,0,
1,0,0,1,1,
1,1,1,1,1)
T_nw <- matrix(c(0,1,1,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0),nrow=n,ncol=n,byrow=TRUE)
act_mat <- calcActivation(T_nw,b,n,K)

# define the parameters for the observation generated from the normal distribution
active_mu <- 0.9
inactive_mu <- 0.5
active_sd <- inactive_sd <- 0.1

# compute the observations matrix
getObsMat(act_mat, active_mu, active_sd, inactive_mu, inactive_sd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{activation }

