\name{calcPrediction}
\alias{calcPredictionLOOCV}
\alias{calcPredictionKfoldCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Predicted Observation.
}
\description{
Calculate the predicted observation of a perturbation experiment. If observations of an experiment are missing this function can be used to determine for a given network the predicted outcome. The missing measurement is predicted from two normal distributions, one for observations coming from active and one coming from inactive genes. The state of the gene is predicted based on the states of its parents.
}
\usage{
calcPredictionLOOCV(obs, delta, b, n ,K, adja, baseline, rem_gene, 
                    rem_k, rem_t=NULL, active_mu, active_sd, inactive_mu, 
                    inactive_sd, mu_type, flag_time_series=FALSE)
calcPredictionKfoldCV(obs, delta, b, n, K, adja, baseline, rem_entries=NULL,
                      rem_entries_vec=NULL, active_mu, active_sd, inactive_mu, 
                      inactive_sd, mu_type, flag_time_series=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
Numeric matrix/array: the observation matrix/array. It can have up to 3 dimensions, where dimension 1 are the network nodes, dimension 2 are the perturbation experiments, and dimension 3 are the time points (if considered).
}
  \item{delta}{
Numeric vector defining the thresholds for each gene to determine its observation to be active or inactive.
}
  \item{b}{
Binary vector representing the perturbation experiments (entry is 0 if gene is inactive in the respective experiment and 1 otherwise). 
}
  \item{n}{
Number of genes in the observation matrix.  
}
  \item{K}{
Number of perturbation experiments.
}
  \item{adja}{
Numeric matrix: the adjacency matrix of the given network.
}
  \item{baseline}{
Vector containing the inferred baseline vectors of each gene.
}
  \item{rem_gene}{
Integer: the index of the gene that is missing.
}
  \item{rem_k}{
Integer: the index of the perturbation experiment that is missing.
}
  \item{rem_t}{
Integer: the index of the time point that is missing.
}
  \item{rem_entries}{
Numeric matrix: each row represents an entry that was removed from the observation matrix, while the 3 columns represent the gene, perturbation experiment and time point, respectively.
}
  \item{rem_entries_vec}{
Numeric vector: contains the entries that were removed in an "absolute form", i.e., if entry (2,1,2) was removed, it will appear in this vector as simply 5.
}
  \item{active_mu}{
Numeric: the average value assumed for observations coming from active nodes. The parameter active_mu and active_sd are used for predicting the observations of the normal distribution of activate genes. This parameter can be either a numeric, a vector, a matrix, or a 3D array, depending on the specified mu_type.
}
  \item{active_sd}{
Numeric: the variation assumed for observations coming from active nodes. The parameter active_mu and active_sd are used for predicting the observations of the normal distribution of activate genes. This parameter can be either a numeric, a vector, a matrix, or a 3D array, depending on the specified mu_type.
}
  \item{inactive_mu}{
Numeric: the average value assumed for observations coming from inactive nodes. The parameter inactive_mu and inactive_sd are used for predicting the observations of the normal distribution of inactivate genes. This parameter can be either a numeric, a vector, a matrix, or a 3D array, depending on the specified mu_type.
}
  \item{inactive_sd}{
  Numeric: the variation assumed for observations coming from inactive nodes. The parameter inactive_mu and inactive_sd are used for predicting the observations of the normal distribution of inactivate genes. This parameter can be either a numeric, a vector, a matrix, or a 3D array, depending on the specified mu_type.
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
  \item{flag_time_series}{
  Logical: specifies whether steady-state (FALSE) or time series data (TRUE) is used.  
}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{loocv}}, \code{\link{kfoldCV}}
}
\examples{
n <- 3 # number of genes
K <- 4 # number of experiments
T_ <- 4 # number of time points

# perturbation vector, entry is 0 if gene is inactivated and 1 otherwise
b <- c(0,1,1, # perturbation exp1: gene 1 perturbed, gene 2,3 unperturbed
       1,0,1, # perturbation exp2: gene 2 perturbed, gene 1,3 unperturbed
       1,1,0, # perturbation exp3....
       1,1,1)

# adjacency matrix       
adja <- matrix(c(0,1,0,
                 0,0,1,
                 0,0,0), nrow=n, ncol=n, byrow=TRUE)

# define node baseline values
baseline <- c(0.75, 0, 0)

# define delta value
delta <- rep(0.75, n)

# define the parameters for the observation generated from the normal distributions
mu_type <- "single"
active_mu <- 0.9
inactive_mu <- 0.5
active_sd <- inactive_sd <- 0.01

#### kfoldCV

# generate random observation matrix
obs <- array(rnorm(n*K*T_), c(n,K,T_))

# define the observationd to be removed, whose values will be predicted
obs[2,4,2] <- NA
obs[3,4,3] <- NA

rem_entries <- which(is.na(obs), arr.ind=TRUE)
rem_entries_vec <- which(is.na(obs))

# compute the predicted observation matrix for the "kfoldCV"
calcPredictionKfoldCV(obs=obs, delta=delta, b=b, n=n, K=K, adja=adja, baseline=baseline, 
                      rem_entries=rem_entries, rem_entries_vec=rem_entries_vec,
                      active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu, 
                      inactive_sd=inactive_sd, mu_type=mu_type, flag_time_series=TRUE)


#### LOOCV
# generate random observation matrix
obs <- matrix(rnorm(n*K), nrow=n, ncol=K)

# define the observationd to be removed
rem_k <- 3
rem_gene <- 2
obs[rem_gene, rem_k] <- NA

# compute the predicted value
calcPredictionLOOCV(obs=obs, delta=delta,b=b, n=n ,K=K, adja=adja, baseline=baseline, 
                    rem_gene=rem_gene, rem_k=rem_k, active_mu=active_mu, active_sd=active_sd,
                    inactive_mu=inactive_mu, inactive_sd=inactive_sd, mu_type=mu_type)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mean squared error }

