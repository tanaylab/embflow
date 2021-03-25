
library("metacell")
scdb_init("scrna_db/", force_reinit=T)
source("scripts/generate_paper_figures/plot_network.r")
#tgconfig::override_params("config/sing_emb.yaml","metacell")


build_sing_emb_wt10_network = function(net_id = "sing_emb_wt10") {
  
  mat_id = "sing_emb_wt10"
  mc_id = "sing_emb_wt10_recolored"
  mgraph_id = "sing_emb_wt10_recolored_logist"
  #net_id = "sing_emb_wt10"
  fig_dir = "figs/sing_emb_wt10.net"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc(mc_id)
  #capacity_var_factor = rep(0.25,ncol(mc@e_gc))
  capacity_var_factor = rep(0.4,ncol(mc@e_gc))
  
  # next define the mc_leak parameter
  
  mc_leak = get_mc_leak_parameter_endo(mc_id = mc_id,leak_emb_endo = 0.12,leak_exe_endo = 0.17)
  
  build_net(mat_id = mat_id,
                    mc_id = mc_id,
                    mgraph_id = mgraph_id,
                    net_id = net_id,
                    fig_dir = fig_dir,
                    age_field = "age_group",mc_leak = mc_leak,
                    capacity_var_factor = capacity_var_factor,
                    k_norm_ext_cost = 1,
                    k_ext_norm_cost = 1,
                    k_ext_ext_cost = 1)
  
  
  
}



build_net = function(mat_id,mc_id,mgraph_id,net_id,fig_dir,
                             age_field = "age_group",
                             mc_leak = NULL,
                             capacity_var_factor = NULL, 
                             t_exp = 1,T_cost = 1e+5,
                             flow_tolerance = 0.01,
                             network_color_ord = NULL,
                             off_capacity_cost1 = 1,
                             off_capacity_cost2 = 1000,
                             k_norm_ext_cost = 2,
                             k_ext_norm_cost = 2,
                             k_ext_ext_cost = 100) {
  
  mat = scdb_mat(mat_id)
  mc  = scdb_mc(mc_id)
  mgraph = scdb_mgraph(mgraph_id)
  md = mat@cell_metadata
  cell_time = md[,age_field]
  names(cell_time) = rownames(md)
  
  if(is.null(mc_leak)) {
    leak = rep(0,max(mc@mc))
  }
  if(is.null(capacity_var_factor)) {
    capacity_var_factor = rep(0.4,max(mc@mc))
  }
  
  if(is.null(mc_leak)) {
    f_extra = mc@colors == "#F6BFCB" | mc@colors == "#7F6874"
    f_endo = mc@colors == "#0F4A9C" | mc@colors == "#EF5A9D" | mc@colors == "#F397C0" | mc@colors == "#c19f70"
    f_pgc = mc@colors == "#FACB12"
    leak = rep(0, max(mc@mc))
    leak[f_extra] = 0.17
    leak[f_endo] = 0.12
  } else {
    leak = mc_leak
  }

  
  mcell_new_mctnetwork(net_id = net_id,
                       mc_id = mc_id,
                       mgraph_id = mgraph_id,
                       cell_time = cell_time)
  mct = scdb_mctnetwork(net_id)

  #computing manifold costs (based on mgraph distances
  mct = mctnetwork_comp_manifold_costs(mct,t_exp=t_exp, T_cost=T_cost)
  message("computed manifold costs")
  
  #generating network structure	
  mct = mctnetwork_gen_network(mct, mc_leak = leak,capacity_var_factor = capacity_var_factor,
                               k_norm_ext_cost = k_norm_ext_cost,k_ext_norm_cost = k_ext_norm_cost,k_ext_ext_cost = k_ext_ext_cost,
                               off_capacity_cost1 = off_capacity_cost1,off_capacity_cost2 = off_capacity_cost2)	
  message("generated network")
  #solving the flow problem
  mct = mctnetwork_gen_mincost_flows(mct, flow_tolerance = flow_tolerance)
  message("solved network flow problem")
  
  #compute propagatation forward and background
  mct = mctnetwork_comp_propagation(mct)
  
  #adding back the object with the network and flows
  scdb_add_mctnetwork(net_id, mct)
  
  mct = scdb_mctnetwork(net_id)
  
  #to plot the "big" network diagram
  if(is.null(network_color_ord)) {
    network_color_ord = mc@color_key$color
  }
  
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mm_mctnetwork_plot_net(mct_id = net_id,
                         fn = sprintf("%s/%s_net.png",fig_dir,net_id),w = 3500, h = 4900,
                         dx_back = 0,
                         colors_ordered=network_color_ord,plot_pdf = F,
                         show_axes = F,
                         show_over_under_flow = F,mc_cex = 1,max_lwd = 15,edge_w_scale = 2e-4,
                         plot_mc_ids = F)


  message("plotted the network")
}



get_mc_leak_parameter_endo = function(mc_id,leak_emb_endo,leak_exe_endo) {
  
  mc = scdb_mc(mc_id)
  
  legc = log2(mc@e_gc + 1e-5)
  
  mc_leak = rep(0,ncol(legc))
  
  # first separation embryonic endoderm (including node/notochord) from meso/-ectoderm 
  x1 = -16
  y1 = -12
  x2 = -12
  y2 = -16
  
  b_emb_endo = (y2 - y1)/(x2 - x1)
  a_emb_endo = (y1*x2 - y2*x1)/(x2 - x1)
  
  f_endo = (legc["Foxa1",]  > a_emb_endo + b_emb_endo*legc["Foxa2",])
  
  mc_leak[f_endo] = leak_emb_endo
  
  # second separation extraembryonic from embryonic endoderm uses Apoe
  x1 = -8.4
  y1 = -14
  x2 = -11
  y2 = -8.4
  
  b_exe_endo = (y2 - y1)/(x2 - x1)
  a_exe_endo = (y1*x2 - y2*x1)/(x2 - x1)
  
  f_exe = (legc["Ttr",]  > a_exe_endo + b_exe_endo*legc["Apoe",])
  
  mc_leak[f_exe] = leak_exe_endo
  
  return(mc_leak)
}







