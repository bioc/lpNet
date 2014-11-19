\name{getObsMat}
\alias{getObsMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Observation Matrix.
}
\description{
The function generates the observation matrix where active/inactive observations are generated from a normal distribution with the average and variation as given in the parameters. This matrix can either be generated from the activation matrix calculated with calcActivation or from the network states caculated with generateTimeSeriesNetStates.
 
}
\usage{
getObsMat(act_mat=NULL, net_states=NULL, active_mu, active_sd, inactive_mu, inactive_sd, mu_type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{act_mat}{
Matrix of 0/1 values called the activation matrix. Rows correspond to genes, columns to experiments. If an entry is 1, it means that the corresponding gene is active in the corresponing experiment and inactive otherwise.
}
  \item{net_states}{
Array of 0/1 values called the network states. Rows correspond to genes, columns to experiments, and the third dimension corresponds to time points. If an entry is 1, it means that the corresponding gene is active in the corresponing experiment and inactive otherwise.
}
  \item{active_mu}{
Numeric: the average value assumed for observations coming from activated nodes. The parameter active_mu and active_sd are used for predicting the observations of the normal distribution of activate genes. This parameter can be either a numeric, a vector, a matrix, or a 3D array, depending on the specified mu_type.
}
  \item{active_sd}{
Numeric: the variation assumed for observations coming from activated nodes. The parameter active_mu and active_sd are used for predicting the observations of the normal distribution of activate genes. This parameter can be either a numeric, a vector, a matrix, or a 3D array, depending on the specified mu_type.
}
  \item{inactive_mu}{
Numeric: the average value assumed for observations coming from inactivated nodes. The parameter inactive_mu and inactive_sd are used for predicting the observations of the normal distribution of inactivate genes. This parameter can be either a numeric, a vector, a matrix, or a 3D array, depending on the specified mu_type.
}
  \item{inactive_sd}{
  Numeric: the variation assumed for observations coming from inactivated nodes. The parameter inactive_mu and inactive_sd are used for predicting the observations of the normal distribution of inactivate genes. This parameter can be either a numeric, a vector, a matrix, or a 3D array, depending on the specified mu_type.
}
  \item{mu_type}{
Character: can have the following values and meanings:
\itemize{
    \item "simple" - the value of active_mu/sd and inactive_mu/sd is independent of the gene/perturbation experiment/time point;
    \item "perGene" - the value of active_mu/sd and inactive_mu/sd depends on the gene;
    \item perGeneExp" - the value of active_mu/sd and inactive_mu/sd depends on the gene and perturbation experiment;
    \item perGeneTime" - the value of active_mu/sd and inactive_mu/sd depends on the gene and time point;
    \item "perGeneExpTime" - the value of active_mu/sd and inactive_mu/sd depends on the gene, perturbation experiment, and time point;
  }
}
}

\value{
Numeric matrix/array: the observation matrix/array. It can have up to 3 dimensions, where dimension 1 are the network nodes, dimension 2 are the perturbation experiments, and dimension 3 are the time points (if considered).
}

\seealso{
\code{\link{calcActivation}}
}
\examples{
n <- 5 # number of genes
K <- 7 # number of knockdowns

# perturbation vector, entry is 0 if gene is inactivated and 1 otherwise
b <- c(0,1,1,1,1, # perturbation exp1: gene 1 perturbed, gene 2-5 unperturbed
       1,0,1,1,1, # perturbation exp2: gene 2 perturbed, gene 1,3,4,5 unperturbed
       1,1,0,1,1, # perturbation exp3....
       1,1,1,0,1,
       1,1,1,1,0,
       1,0,0,1,1,
       1,1,1,1,1)
       
T_nw <- matrix(c(0,1,1,0,0,
                 0,0,0,-1,0,
                 0,0,0,1,0,
                 0,0,0,0,1,
                 0,0,0,0,0), nrow=n, ncol=n, byrow=TRUE)
                 
act_mat <- calcActivation(T_nw, b, n, K)

# define the parameters for the observation generated from the normal distribution
active_mu <- 0.9
inactive_mu <- 0.5
active_sd <- inactive_sd <- 0.1
mu_type <- "single"

# compute the observations matrix
getObsMat(act_mat=act_mat, active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu, inactive_sd=inactive_sd, mu_type=mu_type)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{activation}

