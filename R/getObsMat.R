getObsMat <-
function(act_mat,active_mu,active_sd,inactive_mu,inactive_sd){
  obs_mat <- act_mat
  obs_mat[act_mat==1] <- rnorm(sum(act_mat==1),active_mu,active_sd)
  obs_mat[act_mat==0] <- rnorm(sum(act_mat==0),inactive_mu,inactive_sd)
  return(obs_mat)
}
