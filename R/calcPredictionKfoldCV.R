calcPredictionKfoldCV <-
function(adja,b,n,K,active_mu,active_sd,inactive_mu,inactive_sd){
  # which genes are silenced in removed observation
  act_mat <- calcActivation(adja,b,n,K)
  predict <- getObsMat(act_mat,active_mu,inactive_mu,active_sd,inactive_sd)
  return(predict)
}
