### R code from vignette source 'vignette_lpNet.Rnw'

###################################################
### code chunk number 1: setup1
###################################################
library("lpNet")


###################################################
### code chunk number 2: data1
###################################################
n <- 5 # number of genes
K <- 7 # number of perturbation experiments
obs <- matrix(rnorm(n*K), nrow=n, ncol=K)
delta <- apply(obs, 1, mean)
delta_type <- "perGene"
# perturbation vector (0 if gene inactivated and 1 otherwise)
b <- c(0,1,1,1,1,    # perturbation exp1
       1,0,1,1,1,     # perturbation exp2
       1,1,0,1,1,     # perturbation exp3...
       1,1,1,0,1,
       1,1,1,1,0,
       1,0,0,1,1,
       1,1,1,1,1)


###################################################
### code chunk number 3: lp
###################################################
res1 <- doILP(obs, delta, lambda=1, b, n, K, T_=NULL,
              annot=getEdgeAnnot(n), delta_type)
adja1 <- getAdja(res1, n)


###################################################
### code chunk number 4: lp_pos
###################################################
res2 <- doILP(obs, delta, lambda=1, b, n, K, T_=NULL, 
              annot=getEdgeAnnot(n,allpos=TRUE), delta_type, 
              all.pos=TRUE)
adja2 <- getAdja(res2,n)


###################################################
### code chunk number 5: data_ts
###################################################
T_ <- 4 # number of time points
obs_ts <- array(rnorm(n*K*T_), c(n,K,T_))


###################################################
### code chunk number 6: lp_ts
###################################################
res1 <- doILP(obs_ts, delta, lambda=1, b, n, K, T_, 
              annot=getEdgeAnnot(n), delta_type, 
              flag_time_series=TRUE)
adja1 <- getAdja(res1, n)


###################################################
### code chunk number 7: gaussian
###################################################
active_mu <- 1.1
inactive_mu <- 0.9
active_sd <- inactive_sd <- 0.01
mu_type <- "single"


###################################################
### code chunk number 8: loocv
###################################################
times <- 5 # number of times the removed observation is sampled
annot_node <- seq(1,n) # annotation of the nodes
loocv_res <- loocv(kfold=NULL, times, obs, delta, lambda=1, 
                   b, n, K, T_=NULL, annot=getEdgeAnnot(n), 
                   annot_node, active_mu, active_sd, 
                   inactive_mu, inactive_sd, mu_type,  
                   delta_type)

loocv_res$MSE


###################################################
### code chunk number 9: loocv
###################################################
adja_loocv <- getSampleAdjaMAD(loocv_res$edges_all, n, 
                               annot_node)


###################################################
### code chunk number 10: rangeLambda
###################################################
lambda <- calcRangeLambda(obs, delta, delta_type)
MSE <- Inf
for(lamd in lambda){
  loocv_res <- loocv(kfold=NULL, times, obs, delta, lambda=lamd, 
                     b, n, K, T_=NULL, annot=getEdgeAnnot(n), 
                     annot_node, active_mu, active_sd, 
                     inactive_mu, inactive_sd, mu_type, 
                     delta_type)
  if(loocv_res$MSE<MSE){
		MSE <- loocv_res$MSE
		edges_all <- loocv_res$edges_all
		bestLambda <- lamd
  }
}
adja_bestLambda <- getSampleAdjaMAD(edges_all, n, annot_node)


###################################################
### code chunk number 11: foldcv
###################################################
kfold <- 5
MSE <- Inf
for(lamd in lambda){
  kcv_res <- kfoldCV(kfold, times, obs, delta, lambda=lamd, 
                     b, n, K, T_=NULL, annot=getEdgeAnnot(n),
                     annot_node, active_mu, active_sd, 
                     inactive_mu, inactive_sd, mu_type, 
                     delta_type)
  if(kcv_res$MSE<MSE){
		MSE <- kcv_res$MSE
		edges_all <- kcv_res$edges_all
		bestLambda <- lamd
  }
}
adja_bestLambda <- getSampleAdjaMAD(edges_all, n, annot_node)


###################################################
### code chunk number 12: prior1
###################################################
res3 <- doILP(obs, delta, lambda=1, b, n, K, T_=NULL, 
              annot=getEdgeAnnot(n), delta_type, 
              sourceNode=1, sinkNode=5)

adja3 <- getAdja(res3,n)


###################################################
### code chunk number 13: prior2
###################################################
prior <- list(c("w+_1_2", 1, ">", 1))
res4 <- doILP(obs,delta, lambda=1, b, n, K, T_=NULL, 
              annot=getEdgeAnnot(n), delta_type, 
              prior=prior)
adja4 <- getAdja(res4, n)


###################################################
### code chunk number 14: prior3
###################################################
prior <- list(c("w+_1_2", 1/0.9, ">", 1))
res5 <- doILP(obs, delta, lambda=1, b, n, K, T_=NULL, 
              annot=getEdgeAnnot(n), delta_type, 
              prior=prior)
              
adja5 <- getAdja(res5, n)


###################################################
### code chunk number 15: sahinData
###################################################
data("SahinRNAi2008")
dataStim <- dat.normalized[dat.normalized[ ,17] == 1,-17]
dataUnstim <- dat.normalized[dat.normalized[ ,17] == 0,-17]
# summarize replicates; transpose: rows=genes, cols=experiments
dataSt <- t(summarizeRepl(dataStim, type=mean))
dataUnst <- t(summarizeRepl(dataUnstim, type=mean))


###################################################
### code chunk number 16: sahinData_parameter
###################################################
n <- 16 # number of genes
K <- 16 # number of experiments
annot <- getEdgeAnnot(n)
# perturbation vector;       kd of:
b <- c(0,rep(1,15),          # erbb1
       0,0,rep(1,14),           # erbb1 & 2
       0,rep(1,14),0,	   # erbb1 & 3
       1,0,rep(1,13),0,	 # erbb2 & 3
       rep(1,10),0,1,1,1,1,1,   # IGFR
       rep(1,11),0,1,1,1,1,     # ERalpha
       rep(1,12),0,1,1,1,       # MYC
       rep(1,7),0,rep(1,8),     # AKT1
       rep(1,8),0,rep(1,7),     # MEK1
       rep(1,5),0,rep(1,10),    # CDK2
       rep(1,6),0,rep(1,9),     # CDK4
       rep(1,13),0,1,1,         # CDK6
       1,1,0,rep(1,13),         # p21
       1,1,1,0,rep(1,12),       # p27
       rep(1,4),0,rep(1,11),    # Cyclin D1
       rep(1,14),0,1)           # Cyclin E1


###################################################
### code chunk number 17: sahinData_delta
###################################################
delta <- as.double(dataUnst[ ,1])
delta[11:16] <- mean(dataUnst[ ,1], na.rm=T)


###################################################
### code chunk number 18: sahinData_lp
###################################################
resERBB <- doILP(dataSt[ ,-1], delta, lambda=1.83, 
                 b, n, K, T_=NULL, annot, 
                 delta_type, all.pos=FALSE, 
                 sourceNode=c(1,2,16), sinkNode=10)
adjaERBB <- getAdja(resERBB,n)


