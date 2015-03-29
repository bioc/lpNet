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
loocv(kfold=NULL, times, obs, delta, lambda, b, n, K, T_=NULL, 
      annot, annot_node, active_mu, active_sd, inactive_mu, 
      inactive_sd, mu_type, delta_type, prior=NULL, sourceNode=NULL, 
      sinkNode=NULL, allint=FALSE, allpos=FALSE, flag_ss_v2=FALSE, 
      flag_time_series=FALSE)
kfoldCV(kfold, times, obs, delta, lambda, b, n, K, T_=NULL, 
        annot, annot_node, active_mu, active_sd, inactive_mu, 
        inactive_sd, mu_type, delta_type, prior=NULL,
				sourceNode=NULL, sinkNode=NULL, allint=FALSE, 
        allpos=FALSE, flag_ss_v2=FALSE, flag_time_series=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kfold}{
Integer value of the number "k" in the k-fold cross-calidation.
}
  \item{times}{
Integer: the number of times the cross-validation shall be performed.
}
  \item{obs}{
Numeric matrix/array: the measured observation matrix/array. It can have up to 3 dimensions, where dimension 1 are the network nodes, dimension 2 are the perturbation experiments, and dimension 3 are the time points (if considered).
}
  \item{delta}{
Numeric vector, matrix, or array defining the thresholds to determine an observation active or inactive.
}
  \item{lambda}{
Numeric value defining the penalty parameter lambda. It can range from zero to infinity and it controls the introduction of slack variables in the network inference lp model. 
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
  \item{T_}{
Integer: number of time points in time-series data.
}
  \item{annot}{
Vector of character strings: the annotation of the edges as returned by "getEdgeAnnot".
}
  \item{annot_node}{
Vector of character strings: the annoation of the nodes.
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
  \item{delta_type}{
Character: can have the following values and meanings:
\itemize{
    \item "perGene" - the value of delta depends on the gene;
    \item "perGeneExp" - the value of delta depends on the gene and perturbation experiment;
    \item "perGeneTime" - the value of delta depends on the gene and time point;
    \item "perGeneExpTime" - the value of delta depends on the gene, perturbation experiment, and time point;
  }
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
  \item{flag_time_series}{
  Logical: specifies whether steady-state (FALSE) or time series data (TRUE) is used.  
}
  \item{flag_ss_v2}{
  Logical: if steady-state data is used, specifies whether inactive nodes should contribute to their children activity (FALSE) or not (TRUE).
}
}

\value{
A list of 
 \item{MSE }{The mean squared error (MSE) of predicted and observed measurements of the corresponding cross-validation step.}
 \item{edges_all }{The learned edge weights for each cross-validation step.}
 \item{baseline_all }{The learned baseline weights for each cross-validation step.}
 }

\examples{
n <- 3 # number of genes
K <- 4 # number of experiments
T_ <- 4 # number of time points

annot_node <- seq(1, n)
annot <- getEdgeAnnot(n)

# generate random observation matrix
obs <- array(rnorm(n*K*T_), c(n,K,T_))

baseline <- c(0.75, 0, 0)

# define delta 
delta <- apply(obs, 1, mean, na.rm=TRUE)

# perturbation vector, entry is 0 if gene is inactivated and 1 otherwise
b <- c(0,1,1, # perturbation exp1: gene 1 perturbed, gene 2,3 unperturbed
       1,0,1, # perturbation exp2: gene 2 perturbed, gene 1,3 unperturbed
       1,1,0, # perturbation exp3....
       1,1,1)
       
T_nw <- matrix(c(0,1,0,
                 0,0,1,
                 0,0,0), nrow=n, ncol=n, byrow=TRUE)
colnames(T_nw) <- rownames(T_nw) <- annot_node

## calculate observation matrix with given parameters for the 
# Gaussian distributions for activation and deactivation
active_mu <- 0.95
inactive_mu <- 0.56
active_sd <- inactive_sd <- 0.1

times <- kfold <- 10   # can be increased i.e. to 1000 to produce stable results

mu_type <- "single"
delta_type <- "perGene"

lambda <- 1/10


#### LOOCV
loocv(kfold=NULL, times=times, obs=obs, delta=delta, lambda=lambda, b=b, n=n, K=K, T_=T_, annot=annot,
      annot_node=annot_node, active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu,
      inactive_sd=inactive_sd, mu_type=mu_type, delta_type=delta_type, prior=NULL, sourceNode=NULL,
      sinkNode=NULL, allint=FALSE, allpos=FALSE, flag_time_series=TRUE)

#### K-fold CV
kfoldCV(kfold=kfold, times=times, obs=obs, delta=delta, lambda=lambda, b=b, n=n, K=K, T_=T_, annot=annot,
        annot_node=annot_node, active_mu=active_mu, active_sd=active_sd, inactive_mu=inactive_mu,
        inactive_sd=inactive_sd, mu_type=mu_type, delta_type=delta_type, prior=NULL, sourceNode=NULL,
        sinkNode=NULL, allint=FALSE, allpos=FALSE, flag_time_series=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cross-validation }
