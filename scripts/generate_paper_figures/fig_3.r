source("scripts/generate_paper_figures/plot_network.r")
library("metacell")
library("ggplot2")
library("Matrix")
scdb_init("scrna_db/", force_reinit=T)
scfigs_init("figs")
tgconfig::override_params("config/sing_emb.yaml","metacell")

# wt10 fig3 plots

gen_fig_3_plots = function() {
  
  # The parameters for all plots can be found in the function fig3_parameters
  
  
  if(!dir.exists("figs/paper_figs/fig3")) {
    dir.create("figs/paper_figs/fig3")
  }
  if(!dir.exists("figs/paper_figs/fig_s4")) {
    dir.create("figs/paper_figs/fig_s4")
  }
  
  # The functions below also generate the analogous plots of Figure S4a
  fig3_generate_heatmaps()
  fig3_generate_polygons_cell_type_composition()
  fig3_plot_line_graphs()
  
  fig3_plot_legend_line_graphs()
  
  fig3_plot_color_scale_heatmap()
  
}

fig3_generate_heatmaps = function() {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  mct_id = "sing_emb_wt10"
  mct = scdb_mctnetwork(mct_id)
  
  fig3_param_ls = fig3_parameters()
  
  highlighted_genes = fig3_param_ls$highlighted_genes
  
  
  fig_dir = "figs/paper_figs/fig3"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"age_group"])
  mc_ag_n = mc_ag/rowSums(mc_ag)
  mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
  mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)
  
  late_mcs = which(mc_ag_cn[,13] > 0.2)
  
  df_param = fig3_param_ls$df_param
  
  fig_subdirectory = c("fig3/fig_3a","fig3/fig_3c","fig3/fig_3b","fig_s4/fig_s4a")
  
  for (i in 1:nrow(df_param)) {
    m = df_param$mc[i]
    lfp_thr = df_param$lfp_thr[i]
    min_thr = df_param$min_thr[i]
    
    
    
    fig_dir = sprintf("figs/paper_figs/%s",fig_subdirectory[i])
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    
    p_mc = rep(0,ncol(mc@e_gc))
    p_mc[m] = mc_ag_c[m,ncol(mc_ag_c)]
    
    probs = mctnetwork_propogate_from_t(mct, ncol(mc_ag_c), p_mc)
    mc_prob = t(t(probs$probs)/colSums(probs$probs))
    
    markers = highlighted_genes[[as.character(m)]]
    
    generate_heatmap_along_trajectory(mc_prob = mc_prob,mc = mc,
                                      fig_dir = fig_dir,
                                      tag = paste0("_mc",as.character(m)),
                                      min_max_fold = lfp_thr,
                                      main = sprintf("Metacell %d",m),
                                      min_thr = min_thr,show_rownames = F,
                                      highlighted_genes = markers)
    
    
  }
  
}

fig3_generate_polygons_cell_type_composition = function() {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  mct_id = "sing_emb_wt10"
  mct = scdb_mctnetwork(mct_id)
  
  fig3_param_ls = fig3_parameters()
  
  highlighted_genes = fig3_param_ls$highlighted_genes
  
  mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"age_group"])
  mc_ag_n = mc_ag/rowSums(mc_ag)
  mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
  mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)
  
  late_mcs = which(mc_ag_cn[,13] > 0.2)
  
  df_param = fig3_param_ls$df_param
  
  fig_subdirectory = c("fig3/fig_3a","fig3/fig_3c","fig3/fig_3b","fig_s4/fig_s4a")
  
  for (i in 1:nrow(df_param)) {
    m = df_param$mc[i]
    lfp_thr = df_param$lfp_thr[i]
    min_thr = df_param$min_thr[i]
    
    fig_dir = sprintf("figs/paper_figs/%s",fig_subdirectory[i])
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    
    p_mc = rep(0,ncol(mc@e_gc))
    p_mc[m] = mc_ag_c[m,ncol(mc_ag_c)]
    
    probs = mctnetwork_propogate_from_t(mct, ncol(mc_ag_c), p_mc)
    mc_prob = t(t(probs$probs)/colSums(probs$probs))
    
    markers = highlighted_genes[[as.character(m)]]

    color_key = mc@color_key
    fig_fn = sprintf("%s/polygon_cell_type_composition_traj_mc%d.pdf",fig_dir,m)
    
    
    mc_prob = t(t(mc_prob)/colSums(mc_prob))
    
    ct_prob = tgs_matrix_tapply(t(mc_prob),mc@colors,sum)
    colnames(ct_prob) = c(1:ncol(ct_prob))
    
    col_to_rank = c(1:nrow(color_key))
    names(col_to_rank) = color_key$color
    col_ord = order(col_to_rank[rownames(ct_prob)])
    
    ct_prob = ct_prob[col_ord,]
    
    age_groups = c(1:13)
    
    poly_y = apply(ct_prob,2,cumsum)
    poly_y = rbind(rep(0,ncol(poly_y)),poly_y)
    
    pdf(fig_fn, w=30, h=3,useDingbats = F)
    par(mar = c(0.1,0.01,0.1,0.01))
    plot(NA, xlim=c(0,13), ylim=c(0,1),axes = F,xlab = "",ylab = "")
    for(i in 2:nrow(poly_y)) {
      polygon(x=c(age_groups,rev(age_groups)), 
              y=c(poly_y[i-1,],rev(poly_y[i,])),
              col=rownames(poly_y)[i])
    }
    dev.off()
    
    
    #plot_type_trace(mc_prob,mc,fig_fn = sprintf("%s/%d_polygon.png",fig_dir,m))
    
    
    
  }
  
  
  
}

fig3_plot_line_graphs = function(plot_pdf = T) {
  
  w = 3600
  h = 1400
  
  fig3_param_ls = fig3_parameters()
  
  highlighted_genes = fig3_param_ls$highlighted_genes
  traj_ls = fig3_param_ls$traj_ls
  df_param = fig3_param_ls$df_param
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  mct_id = "sing_emb_wt10"
  mct = scdb_mctnetwork(mct_id)
  
  network_color_ord = mc@color_key$color
  
  mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"age_group"])
  mc_ag_n = mc_ag/rowSums(mc_ag)
  mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
  mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)
  
  fig_subdirectory = c("fig3/fig_3a","fig3/fig_3c","fig3/fig_3b","fig_s4/fig_s4a")
  
  
  for (i in 1:4) {
    
    m = df_param$mc[i]
    lfp_thr = df_param$lfp_thr[i]
    min_thr = df_param$min_thr[i]
    
    fig_dir = sprintf("figs/paper_figs/%s",fig_subdirectory[i])
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    
    p_mc = rep(0,ncol(mc@e_gc))
    p_mc[m] = mc_ag_c[m,ncol(mc_ag_c)]
    
    probs_trans = mctnetwork_propogate_from_t(mct, ncol(mc_ag_c), p_mc)
    mc_prob = t(t(probs_trans$probs)/colSums(probs_trans$probs))
    
    mm_mctnetwork_plot_net(mct_id, fn=sprintf("%s/flow_chart_traj_mc%d.png",fig_dir,m), 
                           propogate = probs_trans$step_m,
                           colors_ordered = network_color_ord,
                           w = w, h =h, mc_cex = 1,
                           edge_w_scale = 1.5*df_param[as.character(m),"edge_w_scale"],
                           fr_scale = 1.5*df_param[as.character(m),"fr_scale"],
                           max_lwd = 1.5*df_param[as.character(m),"max_lwd"],
                           plot_pdf = F)
    
    
    mc_ls = traj_ls[[as.character(m)]]
    
    mc_prob_ls = list()
    probs_trans_ls = list()
    probs_trans_diff = probs_trans
    for(j in 1:length(mc_ls)) {
      
      traj_par = mc_ls[[j]]
      
      types1 = traj_par$types1
      types2 = traj_par$types2
      t_list = traj_par$t1
      fn_tag = traj_par$fn_tag
      
      cmp_probs_trans = mct_extract_subtrajectory(probs_trans_diff,types1,types2,t_list,mct,mc)
      
      probs_trans_diff = cmp_probs_trans$probs_trans_diff
      
      
      mc_prob_ls[[j]] =  cmp_probs_trans$probs_trans_new$probs

    }
    
    markers = highlighted_genes[[as.character(m)]]  
    
    fig_dir = sprintf("%s/genes",fig_dir)
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    
    line_graph_single_gene_along_trajectory(mc_prob_ls = mc_prob_ls,genes = markers,mc = mc,fig_dir = fig_dir,show_sd = F,plot_pdf = plot_pdf)
  }
  
  
}

fig3_plot_legend_line_graphs = function() {
  
  pdf("figs/paper_figs/fig3/legend_line_graphs.pdf",useDingbats = F)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("center", legend =paste0("Traj. ",as.roman(c(1:3))), lty = c("solid","dashed","dotted"))
  dev.off()

  pdf("figs/paper_figs/fig_s4/fig_s4a/legend_line_graphs.pdf",useDingbats = F)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("center", legend =paste0("Traj. ",as.roman(c(1:3))), lty = c("solid","dashed","dotted"))
  dev.off()
}


fig3_plot_color_scale_heatmap = function() {
  
  shades = colorRampPalette(RColorBrewer::brewer.pal(n = 9,name = "PuBu"))(101)
  vals = seq(0,1,length.out = 101)
  cols = shades
  show_vals_ind =  c(1,51,101)
  
  
  pdf(file = "figs/paper_figs/fig3/color_scale_heatmaps.pdf",useDingbats = F)
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')
  
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  dev.off()
  
  pdf(file = "figs/paper_figs/fig_s4/fig_s4a/color_scale_heatmap.pdf",useDingbats = F)
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')
  
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  dev.off()
 
}

fig3_parameters = function() {
  
  
  highlighted_genes = list("1" = c("Utf1","Dnmt3b","Dnmt3a","Eomes","Mixl","T","Tdgf1","Lefty2","Mesp1","Gata6","Hes1","Hand1","Tbx3","Etv2","Kdr","Rspo3","Gata2","Tal1","Runx1",
                                   "Cited2","Lmo2","Car2","Gata1","Cited4","Klf1","Hba-x"),
                           "57" = c("Utf1","Pou3f1","Eomes","Mixl1","Lefty2","Cyp26a1","Mesp1","Tbx3","Msx2","Hand1","Phlda2","Twist1","Crabp1","Gata4","Gata6",
                                    "Myl7","Meis1","Smarcd3","Gata5","Nkx2-5","Myl4","Tnnc1"),
                           "367" = c("Dnmt3b","Tet1","Utf1","Dnmt3a","Epcam","Mycn","Eomes","Fgf5","Foxa2","Chrd","Bmp7","T","Nog","Sp5","Noto","Iqcg","Tppp3","Dnaic1","Ccdc39",
                                     "Cetn4","Fam183b","Wdr69","Foxj1;Rnf157","Rfx3"),
                           "397" = c("Id2","Sox9","Irx3","Otx1","Otx2","Sfrp1","Dnmt3a","Pou3f1","Fgf5","Nefl","Irx2","Ptn","Irx5","Cntfr","Lhx2","Gdf3","Fgf4","Shisa2","Six3",
                                     "Olig3","Oct-3/4;Pou5f1","Tdgf1","Crabp2","Sox1","Pax3","Myc","Utf1","Sox21","Hesx1","Robo3","Dnmt3b","Foxd3","Etv1","Sall3","Fgf15"))
  
  df_param = data.frame(mc = c(1,57,367,397),
                        edge_w_scale = c(5e-5,5e-5,2e-5,5e-5),
                        fr_scale = c(1,1,1,1),
                        max_lwd = c(20,20,20,20),
                        lfp_thr = c(2,2,2,1.5),
                        min_thr = c(-13,-13,-13,-13))
  
  rownames(df_param) = df_param$mc
  
  mc1_ls = list(traj1 = list(types1 = c("#C594BF"),
                             types2 = c("#FBBE92","#c9a997"),
                             t1 = 3,
                             fn_tag = "traj1"),
                traj2 = list(types1 = c("#C594BF"),
                             types2 = c("#FBBE92","#c9a997"),
                             t1 = 4,
                             fn_tag = "traj2"),
                traj3 = list(types1 = c("#C594BF"),
                             types2 = c("#FBBE92","#c9a997"),
                             t1 = 5,
                             fn_tag = "traj3"))
  
  mc57_ls = list(traj1 = list(types1 = c("#C594BF","#DFCDE4"),
                              types2 = c("#53f1fc"),
                              t1 = c(4,5,6),
                              fn_tag = "traj1"),
                 traj2 = list(types1 = c("#C594BF","#DFCDE4"),
                              types2 = c("#53f1fc"),
                              t1 = c(7,8),
                              fn_tag = "traj2"),
                 traj3 = list(types1 = c("#C594BF","#DFCDE4"),
                              types2 = c("#53f1fc"),
                              t1 = c(9,10),
                              fn_tag = "traj3"))
  
  mc367_ls = list(traj1 = list(types1 = c("#c19f70"),
                               types2 = c("#0F4A9C"),
                               t1 = c(3,4,5),
                               fn_tag = "traj1"),
                  traj2 = list(types1 = c("#c19f70"),
                               types2 = c("#0F4A9C"),
                               t1 = c(6),
                               fn_tag = "traj2"),
                  traj3 = list(types1 = c("#c19f70"),
                               types2 = c("#0F4A9C"),
                               t1 = c(7,8),
                               fn_tag = "traj3"))
  
  mc397_ls = list(traj1 = list(types1 = c("#635547"),
                               types2 = c("#65A83E"),
                               t1 = c(4,5,6),
                               fn_tag = "traj1"),
                  traj2 = list(types1 = c("#635547"),
                               types2 = c("#65A83E"),
                               t1 = c(7),
                               fn_tag = "traj2"),
                  traj3 = list(types1 = c("#635547"),
                               types2 = c("#65A83E"),
                               t1 = c(8),
                               fn_tag = "traj3"))
  
  traj_ls = list(mc1 = mc1_ls,mc57 = mc57_ls,mc367 = mc367_ls,mc397 = mc397_ls)
  names(traj_ls) = c(1,57,367,397)
  
  return(list(highlighted_genes = highlighted_genes,traj_ls = traj_ls,df_param = df_param))
}


generate_heatmap_along_trajectory = function(mc_prob,mc,fig_dir,tag,min_max_fold = 2,min_thr = -13,main = "",show_rownames = T,
                                          highlighted_genes = NULL) {
  
  
  legc = log2(mc@e_gc + 1e-5)
  
  mc_prob = t(t(mc_prob)/colSums(mc_prob))
  
  legc_traj = legc %*% mc_prob
  
  min_traj = apply(legc_traj,1,min)
  max_traj = apply(legc_traj,1,max)
  
  f = max_traj- min_traj > min_max_fold & max_traj > min_thr
  
  gnms_f = rownames(legc_traj)[f]
  
  legc_traj_n = (legc_traj[gnms_f,] - min_traj[gnms_f])/(max_traj[gnms_f]-min_traj[gnms_f])
  
  legc_traj_nn = legc_traj_n/rowSums(legc_traj_n)
  
  gene_mean_time = (legc_traj_nn %*% c(1:13))[,1]
  gene_mode = apply(legc_traj_nn,1,which.max)
  gene_rank = 10*gene_mode + gene_mean_time
  
  
  #km_traj = tglkmeans::TGL_kmeans(legc_traj_n,n_gmod,id_column=F)
  #names(km_traj$cluster) = rownames(legc_traj_n)
  
  #cluster_time = tapply(gene_mean_time,km_traj$cluster,mean)
  #cluster_rank = rank(cluster_time)
  
  #gene_rank = 1000*cluster_rank[km_traj$cluster] + rank(gene_mean_time)
  
  gene_ord = order(gene_rank)
  
  legc_traj_n = legc_traj_n[gene_ord,]


  

  
  legc_traj_n = legc_traj_n[c(nrow(legc_traj_n):1),]
  
  #nm_to_pos = c(1:nrow(legc_traj_n))
  nm_to_pos = seq(0,1,length.out = nrow(legc_traj_n))
  names(nm_to_pos) = rownames(legc_traj_n)
  
  
  rownames(legc_traj_n) = substr(rownames(legc_traj_n),1,10)
  colnames(legc_traj_n) = c(1:13)
  
  #gene_clust = km_traj$cluster[gene_ord]
  
  #cluster_gap_position = which(diff(gene_clust) != 0)
  
  shades = colorRampPalette(RColorBrewer::brewer.pal(n = 9,name = "PuBu"))(1000)

  pdf(file = sprintf("%s/heatmap_genes_along_traj%s.pdf",fig_dir,tag),h = 16, w = 10,useDingbats = F)
  par(mar = c(5,12,2,1))
  image(x = t(legc_traj_n),col = shades,axes = F)
  
  ind_even = seq.int(from = 2,by = 2,length.out = (ncol(legc_traj_n) %/% 2))
  ind_odd = seq.int(from = 1,by = 2,length.out = (ncol(legc_traj_n) %/% 2 + ncol(legc_traj_n) %% 2 ))
  axis_pos = seq(0,1,length.out = 13)
  
  axis(side = 1,at = axis_pos[ind_even],labels = c(1:ncol(legc_traj_n))[ind_even],lwd = 0,lwd.ticks = 2,cex.axis = 2,padj = 1)
  axis(side = 1,at = axis_pos[ind_odd],labels = c(1:ncol(legc_traj_n))[ind_odd],lwd = 0,lwd.ticks = 2,cex.axis = 2,padj = 1)
  if(!is.null(highlighted_genes)) {
    highlighted_genes = intersect(names(nm_to_pos),highlighted_genes)
    axis(side = 2,at = nm_to_pos[highlighted_genes],labels = highlighted_genes,lwd = 0,lwd.ticks = 2,cex.axis = 2,las = 2)
  }
  dev.off()

}

mct_extract_subtrajectory = function(probs_trans,types1,types2,t_list,mct,mc) {
  
  n_mc = ncol(mc@e_gc)
  # generate mc_forward and mc_backward matrices from probs_trans 
  # that are used below to propagate p_mc vectors
  
  
  mc_forward = lapply(probs_trans$step_m,function(mat_tr) {
    
    f = rowSums(mat_tr) > 0
    if(sum(f) > 1) {
      mat_tr[f,] = mat_tr[f,]/rowSums(mat_tr[f,])
    } else {
      mat_tr[f,] = mat_tr[f,]/sum(mat_tr[f,])
    }
    
    
    return(mat_tr)
  })
  
  mc_backward = lapply(probs_trans$step_m,function(mat_tr) {
    
    f = colSums(mat_tr) > 0 
    if(sum(f) > 1) {
      mat_tr[,f] = t(t(mat_tr[,f])/colSums(mat_tr[,f]))
    } else {
      mat_tr[,f] = mat_tr[,f]/sum(mat_tr[,f])
    }
    
    
    return(mat_tr)
  })
  
  
  mcs_type1 = which(mc@colors %in% types1)
  mcs_type2 = which(mc@colors %in% types2)
  
  v0 = rep(0,ncol(mc@e_gc))
  
  probs_trans_diff = probs_trans
  probs_trans_new = matrix(0, nrow = nrow(probs_trans$probs), ncol=ncol(probs_trans$probs))
  step_m_new = list()
  for (i in 1:length(probs_trans$step_m)) {
    step_m_new[[i]] = matrix(0,nrow = nrow(probs_trans$probs),ncol = nrow(probs_trans$probs))
  }
  probs_new = matrix(0, nrow = nrow(probs_trans_diff$probs), ncol=ncol(probs_trans_diff$probs))
  
  for (t1 in t_list) {
    transition_t = probs_trans$step_m[[t1]]
    mc_vec2 = v0
    mc_vec2[mcs_type2] = 1
    mc_vec1 = v0
    mc_vec1[mcs_type1] = 1
    f_row = (transition_t %*% mc_vec2)[,1]  > 0 
    f_col = (mc_vec1 %*% transition_t)[1,]  > 0
    mc_list_t1 = intersect(mcs_type1,c(1:n_mc)[f_row])
    mc_list_t2 = intersect(mcs_type2,c(1:n_mc)[f_col])
    
    max_t = ncol(probs_trans_diff$probs)
    step_m = list()
    probs = matrix(0, nrow = nrow(probs_trans_diff$probs), ncol=ncol(probs_trans_diff$probs))
    p_mc_t1 = rep(0,nrow(probs_trans_diff$probs))
    p_mc_t2 = p_mc_t1
    
    mc_vec2 = v0
    mc_vec2[mc_list_t2] = 1
    mc_vec1 = v0
    mc_vec1[mc_list_t1] = 1
    
    p_mc_t1[mc_list_t1] = (probs_trans_diff$step_m[[t1]] %*% mc_vec2)[mc_list_t1,1]
    p_mc_t2[mc_list_t2] = (mc_vec1 %*% probs_trans_diff$step_m[[t1]])[1,mc_list_t2]
    
    probs[,t1] = p_mc_t1
    probs[,t1 + 1] = p_mc_t2
    
    
    
    probs_trans_diff$probs[,t1] = probs_trans_diff$probs[,t1] - p_mc_t1
    probs_trans_diff$probs[,t1 + 1] = probs_trans_diff$probs[,t1 + 1] - p_mc_t2
    
    mat_tr_t1 = Matrix(0,nrow = nrow(probs_trans$probs),ncol = nrow(probs_trans$probs),sparse = T)
    mat_tr_t1[mc_list_t1,mc_list_t2] = probs_trans_diff$step_m[[t1]][mc_list_t1,mc_list_t2]
    step_m[[t1]] = mat_tr_t1
    probs_trans_diff$step_m[[t1]] = probs_trans_diff$step_m[[t1]] - mat_tr_t1
    
    if(t1 > 1) {
      for(i in (t1-1):1) {
        step_m[[i]] = Matrix(t(t(as.matrix(mc_backward[[i]])) * probs[,i+1]), sparse=T)
        
        probs_trans_diff$step_m[[i]] = probs_trans_diff$step_m[[i]] - step_m[[i]]
        
        probs[,i] = as.matrix(mc_backward[[i]]) %*% probs[,i+1]
        
        probs_trans_diff$probs[,i] = probs_trans_diff$probs[,i] - probs[,i]
      }
    }
    if(t1 < max_t - 1) {
      for(i in (t1+2):max_t) {
        step_m[[i-1]] = Matrix(as.matrix(mc_forward[[i-1]]) * probs[,i-1], sparse=T)
        
        probs_trans_diff$step_m[[i-1]] = probs_trans_diff$step_m[[i-1]] - step_m[[i-1]]
        
        probs[,i] = t(probs[,i-1]) %*% as.matrix(mc_forward[[i-1]])
        
        probs_trans_diff$probs[,i] = probs_trans_diff$probs[,i] - probs[,i]
      }
    }
    
    probs_new = probs_new + probs
    
    for (i in 1:length(step_m_new)) {
      step_m_new[[i]] = step_m_new[[i]] + step_m[[i]]
    }
    
  }
  
  probs_trans_new = list(probs=probs_new, step_m=step_m_new)
  
  
  return(list(probs_trans_new = probs_trans_new,probs_trans_diff = probs_trans_diff))
}


line_graph_single_gene_along_trajectory = function(mc_prob_ls,mc,fig_dir,genes,show_sd = T,plot_pdf = T,ylim_max = NULL,ylim_min = NULL,add_gridline = F,col = "black",tag = NULL) {
  
  legc = log2(mc@e_gc + 1e-5)
  
  out_ls = lapply(mc_prob_ls,function(mc_prob) {
    mc_prob = t(t(mc_prob)/colSums(mc_prob))
    
    
    
    feat_traj = comp_e_gc_along_trajectory(mc_prob,legc)
    
    e_gc_mean = feat_traj$e_gc_mean
    e_gc_sd = feat_traj$e_gc_sd
    
    
    df_mean = as.data.frame(as.table(e_gc_mean))
    colnames(df_mean) = c("gene","age_group","e_gc")
    df_sd = as.data.frame(as.table(e_gc_sd))
    colnames(df_sd) = c("gene","age_group","sd_e_gc")
    
    if(!identical(df_mean[,1:2],df_sd[,1:2])) {
      stop("dataframes df_mean and df_all are not equal")
    }
    df_all = merge(df_mean,df_sd,by = c("gene","age_group"))
    
    df_all$age_group = as.numeric(df_all$age_group)
    df_all$y_min = df_all$e_gc - df_all$sd_e_gc
    df_all$y_max = df_all$e_gc + df_all$sd_e_gc
    
    if(!dir.exists(fig_dir)) {
      dir.create(fig_dir)
    }
    
    if(is.null(ylim_max)) {
      ylim_max = max(feat_traj$e_gc_mean)
    }
    
    if(is.null(ylim_min)) {
      ylim_min = min(feat_traj$e_gc_mean)
    }
    return(list(df_all = df_all,ylim_max = ylim_max,ylim_min = ylim_min))
  })
  
  df_all_ls = lapply(c(1:length(out_ls)),function(i) {
    a = out_ls[[i]]
    df_all = a$df_all
    
    traj_n = paste0(rep("I",i),collapse = "")
    
    df_all$traj = as.character(rep(traj_n,nrow(df_all)))
    
    return(df_all)
  })
  df_all = do.call(rbind,df_all_ls)
  
  ylim_min_ls = lapply(out_ls,function(a) {
    return(a$ylim_min)
  })
  ylim_min = min(unlist(ylim_min_ls))
  
  ylim_max_ls = lapply(out_ls,function(a) {
    return(a$ylim_max)
  })
  ylim_max = max(unlist(ylim_max_ls))
  
  for (gene in genes) {
    df_f = df_all[df_all$gene == gene,]
    gene_name = gsub(";","_",gene)
    gene_name = gsub("/","_",gene_name)
    if(!show_sd) {
      out = ggplot(data = df_f,aes(x = age_group,y = e_gc,group = traj)) +
        xlab("") +
        ylab("") +
        scale_linetype_manual(values = c("solid","dashed","dotted","twodash")) +
        geom_line(size = 1,aes(linetype = traj)) +
        geom_point(size = 2) + 
        ggtitle(label = gene) + 
        theme(legend.position="none") +
        theme(panel.background = element_rect(fill = "white",color = "black")) +
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              plot.title = element_text(size = 60,hjust = 0.5))
      
    } else {
      out = ggplot(data = df_f,aes(x = age_group,y = e_gc,group = traj)) +
        xlab("") +
        ylab("") +
        geom_line(size = 1,aes(group = traj)) + 
        geom_errorbar(aes(ymax = y_max,ymin = y_min),size = 1.3,width = 0.5) +
        ggtitle(label = gene) + theme(legend.position="none") + 
        theme(panel.background = element_rect(fill = "white",color = "black")) +
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              plot.title = element_text(size = 100,hjust = 0.5))
      
    }
    
    if(add_gridline) {
      out = out +
        geom_vline(xintercept = 5, linetype = "dashed",size = 2) +
        geom_vline(xintercept = 10, linetype = "dashed",size = 2)
    }
    
    if(is.null(tag)) {
      fn = sprintf("%s/%s",fig_dir,gene_name)
    } else {
      fn = sprintf("%s/%s_%s",fig_dir,gene_name,tag)
    }
    
    if(show_sd){
      fn =  sprintf("%s_with_sd.pdf",fn)
    } else {
      fn =  sprintf("%s_without_sd.pdf",fn)
    }
    

    
    if (plot_pdf) {
      ggsave(plot = out,filename = fn,width = 10,height = 7)
    } else {
      fn = gsub(pattern = ".pdf",replacement = ".png",x = fn)
      ggsave(plot = out,filename = fn,width = 10,height = 7,bg = "transparent")
    }
    
  }
  
  
}


comp_e_gc_along_trajectory = function(mc_prob,e_gc) {
  
  mc_prob = t(t(mc_prob)/colSums(mc_prob))
  
  if (is.vector(e_gc)) {
    e_gc_mean = e_gc %*% mc_prob
    
    e_gc_sd = apply(X = mc_prob,MARGIN = 2,FUN = function(p) {
      
      e_gc_m = sum(e_gc * p) 
      var_e_gc = e_gc - e_gc_m
      var_e_gc = var_e_gc ^ 2
      var_e_gc = sum(var_e_gc * p)
      sd_e_gc = sqrt(var_e_gc)
      
      return(sd_e_gc)
    })
  } else {
    
    e_gc_mean = e_gc %*% mc_prob
    
    e_gc_sd = apply(X = mc_prob,MARGIN = 2,FUN = function(p) {
      e_gc_m = (e_gc %*% p)[,1]
      var_e_gc = e_gc - e_gc_m
      var_e_gc = var_e_gc ^ 2
      var_e_gc = (var_e_gc %*% p)[,1]
      sd_e_gc = sqrt(var_e_gc)
      
      return(sd_e_gc)
    })
  }
  
  return(list(e_gc_mean = e_gc_mean, e_gc_sd = e_gc_sd))
}


