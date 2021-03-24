
mctnetwork_get_egc_on_cluster_transition = function(mct, min_time, max_time, type1, type2, mc_type=NULL) {
  mc = scdb_mc(mct@mc_id)
  e_gc = mc@e_gc
  net = mct@network
  
  if(is.null(mc_type)) {
    mc_type = mc@colors
    names(mc_type) = as.character(1:length(mc_type))
  }
  
  #	flow_mm = mctnetwork_get_flow_mat(mct, time, max_time=time)
  
  f_t = net$time1 >= min_time & net$time2 <= max_time &
    net$time1 == net$time2-1 &
    net$type1 != "growth" & net$type2!="growth"
  
  net = net[f_t,]
  f_types = mc_type[as.numeric(net$mc1)]==type1 & mc_type[as.numeric(net$mc2)]==type2
  
  net = net[f_types,]
  
  src_mc_wgt = tapply(net$flow, net$mc1, sum)	
  targ_mc_wgt = tapply(net$flow, net$mc2, sum)
  src_mc_wgt_n = as.vector(src_mc_wgt/sum(src_mc_wgt))
  names(src_mc_wgt_n) = names(src_mc_wgt)
  targ_mc_wgt_n = as.vector(targ_mc_wgt/sum(targ_mc_wgt))
  names(targ_mc_wgt_n) = names(targ_mc_wgt)
  
  src_e_gc = colSums(t(e_gc[,names(src_mc_wgt_n)]) * src_mc_wgt_n)
  targ_e_gc = colSums(t(e_gc[,names(targ_mc_wgt_n)]) * targ_mc_wgt_n)
  
  return(data.frame(src = src_e_gc, targ = targ_e_gc, lf = log2(1e-5+targ_e_gc)-log2(1e-5+src_e_gc)))
}
