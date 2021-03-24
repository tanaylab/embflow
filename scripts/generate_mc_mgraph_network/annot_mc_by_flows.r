
cluster_metacells_by_flow = function(mct_id = "sing_emb_wt10",K = 65) {
  
  clst_flows = mctnetwork_clust_flows(mct_id, K)
  
  fclst = clst_flows$clust
  mc_ord = clst_flows$hc$order
  
  mc_cluster = data.frame(mc = mc_ord,cluster = fclst[mc_ord],mc_rank = c(1:length(mc_ord)))
  return(mc_cluster) 
}