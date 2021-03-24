library("metacell")
tgconfig::override_params("config/sing_emb.yaml","metacell")
source("scripts/generate_mc_mgraph_network/gen_network.r")
scdb_init("scrna_db/stability_analysis/",force_reinit = T)

gen_parameter_stability_analysis = function() {
  # first generate networks for all the parameters

  if(0) {
    generate_mgraph_and_network_for_param_values(1)
    generate_mgraph_and_network_for_param_values(2)
    generate_mgraph_and_network_for_param_values(3)
  }

  # generate a summary list of the flows and plot them
  gen_net_id_ls()
  
  # next generate color summary
  gen_col_summary_all()
}

generate_mgraph_and_network_for_param_values = function(n_split = 1) {
  fig_dir = "figs/paper_figs/fig_s_parameter_stability"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  
  
  mat_id = "sing_emb_wt10"
  mc_id = "sing_emb_wt10_recolored"
  gset = scdb_gset("sing_emb_wt10")
  feat_genes = names(gset@gene_set)
  mgraph_id_orig = "sing_emb_wt10_recolored_logist"
  mc = scdb_mc(mc_id)
  mc_leak = get_mc_leak_parameter_endo(mc_id = mc_id,leak_emb_endo = 0.12,leak_exe_endo = 0.17)
  
  a = 1
  # parameters of interest assigned to their current standard value
  
  logist_loc = 1
  logist_scale = 0.2
  logist_eps = 4e-5
  max_d_fold = 3
  t_exp = 1
  cap_var_factor = 0.4
  capacity_var_factor = rep(cap_var_factor,ncol(mc@e_gc))
  off_capacity_cost1 = 1
  off_capacity_cost2 = 1000
  max_degree_mc = 4
  tgconfig::set_param("mcell_mc2d_max_confu_deg",max_degree_mc,"metacell")

  if(n_split == 1) {
    # start with t_exp
    # sample ten values
    param_ls = 2^(c((-3:6)))
    
    for(i in 1:length(param_ls)) {
      print(i)
      t_exp = param_ls[i]
      net_id = sprintf("t_exp_%d",i)
      
      build_net(mat_id = mat_id,
                        mc_id = mc_id,
                        mgraph_id = mgraph_id_orig,
                        net_id = net_id,
                        fig_dir = fig_dir,
                        age_field = "age_group",mc_leak = mc_leak,
                        capacity_var_factor = capacity_var_factor,
                        k_norm_ext_cost = 1,
                        k_ext_norm_cost = 1,
                        k_ext_ext_cost = 1,
                        t_exp = t_exp)
    }
    t_exp = 1
    
    # next max degree mc
    param_ls = c(3,4,6,8,10,15,20,50)
    
    for(i in 1:length(param_ls)) {
      
      print(i)
      max_degree_mc = param_ls[i]
      tgconfig::set_param("mcell_mc2d_max_confu_deg",max_degree_mc,"metacell")
      mgraph_id = sprintf("max_degree_mc_%d",i)
      net_id = mgraph_id
      
      mgraph = mc2d_comp_mgraph_param(mc, feat_genes, logist_loc, logist_scale, logist_eps, max_d_fold)
      scdb_add_mgraph(id = mgraph_id,mgraph = tgMCManifGraph(mc_id = mc_id,mgraph = mgraph))
      
      build_singemb_net(mat_id = mat_id,
                        mc_id = mc_id,
                        mgraph_id = mgraph_id,
                        net_id = net_id,
                        fig_dir = fig_dir,
                        age_field = "age_group",mc_leak = mc_leak,
                        capacity_var_factor = capacity_var_factor,
                        k_norm_ext_cost = 1,
                        k_ext_norm_cost = 1,
                        k_ext_ext_cost = 1,
                        t_exp = t_exp)
      
    }
    tgconfig::set_param("mcell_mc2d_max_confu_deg",4,"metacell")
    
    
  }
 
  if (n_split == 2) {
    # next logistic loc 
    param_ls = c(0,0.25,0.5,0.75,1,1.5,2,2.5)
    
    for(i in 1:length(param_ls)) {
      
      print(i)
      logist_loc = param_ls[i]
      mgraph_id = sprintf("logist_loc_%d",i)
      net_id = mgraph_id
      
      mgraph = mc2d_comp_mgraph_param(mc, feat_genes, logist_loc, logist_scale, logist_eps, max_d_fold)
      scdb_add_mgraph(id = mgraph_id,mgraph = tgMCManifGraph(mc_id = mc_id,mgraph = mgraph))
      
      build_singemb_net(mat_id = mat_id,
                        mc_id = mc_id,
                        mgraph_id = mgraph_id,
                        net_id = net_id,
                        fig_dir = fig_dir,
                        age_field = "age_group",mc_leak = mc_leak,
                        capacity_var_factor = capacity_var_factor,
                        k_norm_ext_cost = 1,
                        k_ext_norm_cost = 1,
                        k_ext_ext_cost = 1,
                        t_exp = t_exp)
      
    }
    logist_loc = 1
    
    # next logistic scale
    param_ls = c(0.10,0.15,0.20,0.25,0.30,0.35,0.45,0.60,0.80)
    
    for(i in 1:length(param_ls)) {
      
      print(i)
      logist_scale = param_ls[i]
      mgraph_id = sprintf("logist_scale_%d",i)
      net_id = mgraph_id
      
      mgraph = mc2d_comp_mgraph_param(mc, feat_genes, logist_loc, logist_scale, logist_eps, max_d_fold)
      scdb_add_mgraph(id = mgraph_id,mgraph = tgMCManifGraph(mc_id = mc_id,mgraph = mgraph))
      
      build_singemb_net(mat_id = mat_id,
                        mc_id = mc_id,
                        mgraph_id = mgraph_id,
                        net_id = net_id,
                        fig_dir = fig_dir,
                        age_field = "age_group",mc_leak = mc_leak,
                        capacity_var_factor = capacity_var_factor,
                        k_norm_ext_cost = 1,
                        k_ext_norm_cost = 1,
                        k_ext_ext_cost = 1,
                        t_exp = t_exp)
      
    }
    logist_scale = 0.2
  }
  
 
  if(n_split == 3) {
    
    # next follows capacity variance factor
    
    param_ls = c(0.01,0.05,0.1,0.25,0.4,0.6,0.8)
    
    for(i in 1:length(param_ls)) {
      
      print(i)
      cap_var_factor = param_ls[i]
      capacity_var_factor = rep(cap_var_factor,ncol(mc@e_gc))
      net_id = sprintf("cap_var_factor_%d",i)
      
      build_singemb_net(mat_id = mat_id,
                        mc_id = mc_id,
                        mgraph_id = mgraph_id_orig,
                        net_id = net_id,
                        fig_dir = fig_dir,
                        age_field = "age_group",mc_leak = mc_leak,
                        capacity_var_factor = capacity_var_factor,
                        k_norm_ext_cost = 1,
                        k_ext_norm_cost = 1,
                        k_ext_ext_cost = 1,
                        t_exp = t_exp)
      
    }
    capacity_var_factor = rep(0.25,ncol(mc@e_gc))
    
    
    # next follows off_capacity_cost1
    
    param_ls = c(0,1,2,5,10,20,50,100)
    
    for(i in 1:length(param_ls)) {
      
      print(i)
      off_capacity_cost1 = param_ls[i]
      net_id = sprintf("off_capacity_cost1_%d",i)
      
      build_singemb_net(mat_id = mat_id,
                        mc_id = mc_id,
                        mgraph_id = mgraph_id_orig,
                        net_id = net_id,
                        fig_dir = fig_dir,
                        age_field = "age_group",mc_leak = mc_leak,
                        capacity_var_factor = capacity_var_factor,
                        k_norm_ext_cost = 1,
                        k_ext_norm_cost = 1,
                        k_ext_ext_cost = 1,
                        t_exp = t_exp,
                        off_capacity_cost1 = off_capacity_cost1,
                        off_capacity_cost2 = off_capacity_cost2)
      
    }
    off_capacity_cost1 = 1
    
    # next follows off_capacity_cost2
    
    param_ls = c(1,10,20,50,100,200,500,1000,2000,5000)
    
    for(i in 1:length(param_ls)) {
      
      print(i)
      off_capacity_cost2 = param_ls[i]
      net_id = sprintf("off_capacity_cost2_%d",i)
      
      build_singemb_net(mat_id = mat_id,
                        mc_id = mc_id,
                        mgraph_id = mgraph_id_orig,
                        net_id = net_id,
                        fig_dir = fig_dir,
                        age_field = "age_group",mc_leak = mc_leak,
                        capacity_var_factor = capacity_var_factor,
                        k_norm_ext_cost = 1,
                        k_ext_norm_cost = 1,
                        k_ext_ext_cost = 1,
                        t_exp = t_exp,
                        off_capacity_cost1 = off_capacity_cost1,
                        off_capacity_cost2 = off_capacity_cost2)
      
    }
    off_capacity_cost2 = 1000
    
  }
  
  
 
}


gen_net_id_ls = function() {
  
  scdb_init("scrna_db/",force_reinit = T)
  mc = scdb_mc("sing_emb_wt10_recolored")
  network_color_ord = mc@color_key$color
  scdb_init("scrna_db/stability_analysis/",force_reinit = T)
  
  net_id_ls = c()
  all_param_ls = list()
  # start with t_exp
  # sample ten values
  param_ls = 2^(c((-3:6)))
  net_ids = paste0("t_exp_",c(1:length(param_ls)))
  net_id_ls = c(net_id_ls,net_ids)
  all_param_ls[[1]] = param_ls 
  
  # next logistic loc 
  param_ls = c(0,0.25,0.5,0.75,1,1.5,2,2.5)
  net_ids = paste0("logist_loc_",c(1:length(param_ls)))
  net_id_ls = c(net_id_ls,net_ids)
  all_param_ls[[2]] = param_ls
  
  # next logistic scale
  param_ls = c(0.10,0.15,0.20,0.25,0.30,0.35,0.45,0.60,0.80)
  net_ids = paste0("logist_scale_",c(1:length(param_ls)))
  net_id_ls = c(net_id_ls,net_ids)
  all_param_ls[[3]] = param_ls
  
  # next max degree mc
  param_ls = c(3,4,6,8,10,15,20,50)
  net_ids = paste0("max_degree_mc_",c(1:length(param_ls)))
  net_id_ls = c(net_id_ls,net_ids)
  all_param_ls[[4]] = param_ls
  
  # next follows capacity variance factor
  param_ls = c(0.01,0.05,0.1,0.25,0.4,0.6,0.8)
  net_ids = paste0("cap_var_factor_",c(1:length(param_ls)))
  net_id_ls = c(net_id_ls,net_ids)
  all_param_ls[[5]] = param_ls
  
  # next follows off_capacity_cost1
  param_ls = c(0,1,2,5,10,20,50,100)
  net_ids = paste0("off_capacity_cost1_",c(1:length(param_ls)))
  net_id_ls = c(net_id_ls,net_ids)
  all_param_ls[[6]] = param_ls
  
  # next follows off_capacity_cost2
  param_ls = c(1,10,20,50,100,200,500,1000,2000,5000)
  net_ids = paste0("off_capacity_cost2_",c(1:length(param_ls)))
  net_id_ls = c(net_id_ls,net_ids)
  all_param_ls[[7]] = param_ls
  
  
  
  ignore_ls = c()
  ignore_flow = c()
  for(net_id in net_id_ls) {
    print(net_id)
    mct = scdb_mctnetwork(net_id)
    
    tot_flow = sum(mct@network$flow)
    
    if((tot_flow < 20) | (tot_flow > 30)) {
      ignore_ls = c(ignore_ls,net_id)
      ignore_flow = c(ignore_flow,tot_flow)
    }
    
  }
  names(ignore_flow) = ignore_ls
  net_id_ls_f = setdiff(net_id_ls,ignore_ls)

  param_nms = c("t_exp","logist_loc","logist_scale","max_degree_mc","cap_var_factor","off_capacity_cost1","off_capacity_cost2")
  
  for (i in 1:length(all_param_ls)) {
    
    n_param = length(all_param_ls[[i]])
    
    param_nm = param_nms[i]
    
    net_ids = paste(param_nm,c(1:n_param),sep = "_")
    
    f = net_ids %in% net_id_ls_f
    
    all_param_ls[[i]] = all_param_ls[[i]][f]
  }
  
  names(all_param_ls) = param_nms
  
  ls_of_param = all_param_ls
  
  save(net_id_ls,file = "data/parameter_stability_analysis/net_id_ls.Rda")
  save(net_id_ls_f,file = "data/parameter_stability_analysis/net_id_ls_f.Rda")
  save(ls_of_param,file = "data/parameter_stability_analysis/ls_of_param.Rda")
}

gen_col_summary_all = function() {
  
  param_ls = c("t_exp","logist_loc","logist_scale","max_degree_mc","cap_var_factor","off_capacity_cost1","off_capacity_cost2")
  
  load(file = "data/parameter_stability_analysis/net_id_ls.Rda")
  load("data/parameter_stability_analysis/net_id_ls_f.Rda")
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  network_color_ord = mc@color_key$color
  
  for (param in param_ls) {
    print(param)
    net_ids = grep(pattern = param,x = net_id_ls_f,v = T)
    
    net_id = net_ids[1]
    #load(sprintf("data/parameter_stability_analysis/%s.Rda",net_id))
    col_trans_per_t = param_gen_col_summary(net_id)
    
    
    col_dist_all = lapply(col_trans_per_t,function(col_dist) {
      return(list(col_dist))
    })
    
    for (i in 2:length(net_ids)) {
      print(i)
      net_id = net_ids[i]
      #load(sprintf("data/parameter_stability_analysis/%s.Rda",net_id))
      col_trans_per_t = param_gen_col_summary(net_id)
      print(i)
      for (t in 1:12) {
        col_dist_ls_t = col_dist_all[[t]]
        col_dist_ls_t = c(col_dist_ls_t,list(col_trans_per_t[[t]]))
        col_dist_all[[t]] = col_dist_ls_t
      }
      
    }
    save(col_dist_all,file = sprintf("data/parameter_stability_analysis/%s_all.Rda",param))
    
  }
  
}

param_gen_col_summary = function(net_id) {
  scdb_init("scrna_db/stability_analysis/",force_reinit = T)
  mct = scdb_mctnetwork(net_id)
  scdb_init("scrna_db/",force_reinit = T)
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  all_colors = unique(mc@colors)
  
  cls_f = names(mc@mc)
  
  cls_spl = split(cls_f,mat@cell_metadata[cls_f,"age_group"])
  
  col_dist_zero = matrix(0,nrow = ncol(mc@e_gc),ncol = length(all_colors))
  rownames(col_dist_zero) = c(1:ncol(mc@e_gc))
  colnames(col_dist_zero) = all_colors
  
  col_dist_per_t = lapply(cls_spl,function(cls) {
    
    col_dist = col_dist_zero
    col_dist_tmp = table(mc@mc[cls],mc@colors[mc@mc[cls]])
    col_dist_tmp = col_dist_tmp/rowSums(col_dist_tmp)
    
    col_dist[rownames(col_dist_tmp),colnames(col_dist_tmp)] = col_dist_tmp
    
    return(col_dist)
  })
  
  
  col_trans_per_t = lapply(c(1:length(mct@mc_forward)),function(i) {
    
    mc_trans = mct@mc_t_infer[,i]*mct@mc_forward[[i]]
    
    col_trans = t(col_dist_per_t[[i]]) %*% mc_trans %*% col_dist_per_t[[i+1]]
    return(col_trans)
  })
  
  #save(col_trans_per_t,file = sprintf("data/parameter_stability_analysis/%s.Rda",net_id))
  return(col_trans_per_t)
}
