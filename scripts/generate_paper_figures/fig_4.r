library("shape")
library("metacell")
scdb_init("scrna_db/", force_reinit=T)
scfigs_init("figs")
tgconfig::override_params("config/sing_emb.yaml","metacell")

source("scripts/generate_paper_figures/plot_vein.r")

gen_fig_4_plots = function() {
  
  if(!dir.exists("figs/paper_figs/fig4")) {
    dir.create("figs/paper_figs/fig4")
  }
  # epiblast generate heatmaps
  fig_4a()
  fig_4b()
  fig_4c()
  fig_4d()
  fig_4e()
  fig_4f()
  fig_4g()
  fig_4h()
  fig_4i()
  
}

fig_4a  = function() {
  
  if(!dir.exists("figs/paper_figs/fig4")) {
    dir.create("figs/paper_figs/fig4")
  }
  if(!dir.exists("figs/paper_figs/fig4/fig_4a")) {
    dir.create("figs/paper_figs/fig4/fig_4a")
  }
  plot_focals(fig_dir = "figs/paper_figs/fig4/fig_4a",focals = "Epiblast",plot_pdf = T)
  
}


fig_4b = function(plot_pdf = T) {
  
  # calculates and plots epiblast direct cell type fates according to the flows for each time point
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  mct = scdb_mctnetwork("sing_emb_wt10")
  
  epi_mcs = which(mc@colors == "#635547")
  
  mc_list = epi_mcs
  fn = "figs/paper_figs/fig4/fig_4b.pdf"
  time_points = c(1:8)
  
  # mc_t_infer is a matrix of dimension #metacells times #time_points
  # each column is the infered time distribution over metacells from the flows for that time point
  # below we normalize the subset of epiblast metacells
  
  mc_ag_n_f = t(t(mct@mc_t_infer[mc_list,time_points])/colSums(mct@mc_t_infer[mc_list,time_points]))
  rownames(mc_ag_n_f) = mc_list
  colnames(mc_ag_n_f) = time_points
  
  # mct@mc_forward is a list of 12 matrices (number of time points - 1) with nrow = ncol = number of metacells
  # each matrix has rowSums equal to 1 or (0) and gives for each metacell the normalized forward outgoing flow for this time point
  t = time_points[1]
  mc_fate = (mc_ag_n_f[as.character(mc_list),as.character(t)] %*% mct@mc_forward[[t]][mc_list,])[1,]
  
  for (i in time_points[-1]) {
    
    fate_tmp = (mc_ag_n_f[as.character(mc_list),as.character(i)] %*% mct@mc_forward[[i]][mc_list,])[1,]
    mc_fate = cbind(mc_fate,fate_tmp)
  }
  
  fate_ct = tgs_matrix_tapply(t(mc_fate),mc@colors,sum)
  colnames(fate_ct) = time_points
  fate_ct = fate_ct[rowSums(fate_ct)> 0,]
  if(plot_pdf) {
    pdf(fn,useDingbats = F)
  } else {
    png(fn)
  }
  barplot(fate_ct,col = rownames(fate_ct),yaxt = 'n')
  dev.off()
  
  
}


fig_4c = function() {
  
  
  lfp_thr = 2
  min_thr = -14
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  mat_n = t(t(mat@mat)/colSums(mat@mat))
  #mct = scdb_mctnetwork("sing_emb_wt10_cap040")
  fig_dir = "figs/paper_figs/fig4/fig_4c"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  epi_mcs = which(mc@colors == "#635547")
  mc_list = epi_mcs
  time_points = c(1:8)
  
  # compute bulk time dynamics in nascent mesoderm
  bad_genes =  read.table(file = "data/external_data/sing_emb_wt10.bad_genes.txt",sep ="\t",stringsAsFactors = F)
  bad_genes = bad_genes$x
  
  cls = names(mc@mc)[(mc@mc %in% mc_list) & (mat@cell_metadata[names(mc@mc),"age_group"] %in% time_points)]
  
  egc_t = t(tgs_matrix_tapply(mat_n[,cls],mat@cell_metadata[cls,"age_group"],mean))
  rownames(egc_t) = rownames(mat_n)
  
  legc = log2(egc_t + 1e-5)
  
  min_legc = apply(legc,1,min)
  max_legc = apply(legc,1,max)
  
  f = max_legc > min_thr & (max_legc - min_legc) > lfp_thr
  variable_genes = rownames(legc)[f]
  
  epiblast_genes = rownames(mc@mc_fp)[rowMeans(log2(mc@mc_fp[,epi_mcs])) > 0.5]
  
  constitutive_genes = setdiff(epiblast_genes,c(bad_genes,variable_genes))
  
  epiblast_plot_bulk_dynamics(mc_list = epi_mcs,time_points = time_points,fn = "figs/paper_figs/fig4/fig_4c/epiblast_bulk_dynamics_constitutive_genes.pdf",
                              feat_genes = constitutive_genes,cluster_genes = T,n_clust = 1,highlighted_genes = constitutive_genes,
                              w = 12,h = 18)
  
  
  highlighted_genes = c("Dppa4","Tdgf1","Smad7","Apln","Irx3","Foxd3")
  epiblast_plot_bulk_dynamics(mc_list = epi_mcs,time_points = time_points,fn = "figs/paper_figs/fig4/fig_4c/epiblast_bulk_dynamics_variable_genes.pdf",
                              feat_genes = variable_genes,cluster_genes = T,n_clust = 5,highlighted_genes = variable_genes,
                              w = 12,h = 18)
  
}

fig_4d = function() {
  
  
  mct = scdb_mctnetwork("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  confu = mctnetwork_get_flow_mat(mct, -2)
  diag(confu) = 0
  
  mat = scdb_mat("sing_emb_wt10")
  mc_ag = table(mc@mc, mat@cell_metadata[names(mc@mc), "age_group"])
  mc_t = apply(mc_ag,1, function(x) sum((1:13)*x)/sum(x))
  
  genes = c("Eomes","Tdgf1","Irx5")
  
  incl_colors = c("#635547","#DABE99","#9e6762","#65A83E","#647a4f","#354E23")
  
  for(gene in genes) {
    plot_time_gene_mc_flow(confu = confu,mc = mc,mc_t = mc_t,gene = gene,incl_colors = incl_colors,fig_dir = "figs/paper_figs/fig4/fig_4d",max_t = 13,fig_pref = "epiblast",plot_pdf = T)
  }
  
  
}


fig_4e = function(plot_pdf = T) {
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat_n = t(t(mat@mat)/colSums(mat@mat))
  
  fig_dir = "figs/paper_figs/fig4/fig_4e"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  # first endoderm
  genes = c("Lefty1","Cer1")
  
  included_colors = c("#F6BFCB","#F397C0","#EF5A9D")
  
  mc_f = which(mc@colors %in% included_colors)
  cls_f = names(mc@mc)[mc@mc %in% mc_f]
  
  for (gene in genes) {
    legc_bulk = log2(tapply(mat_n[gene,cls_f],mat@cell_metadata[cls_f,"age_group"],mean) + 1e-5)
    
    if(plot_pdf) {
      pdf(sprintf("%s/endoderm_bulk_expr_%s.pdf",fig_dir,gene),width = 5, h = 3.5,useDingbats = F)
      plot(x = c(1:13),y = legc_bulk,xlab = "Age group",ylab = "",main = gene,cex.main = 2,pch = 19,type = "l",lwd = 4,col = "black",ylim = c(-17,-9))
      dev.off()
    } else {
      png(sprintf("%s/endoderm_bulk_expr_%s.png",fig_dir,gene),width = 800, h = 560)
      plot(x = c(1:13),y = legc_bulk,xlab = "Age group",ylab = "",main = gene,cex.main = 2,pch = 19,type = "l",lwd = 4,col = "black",ylim = c(-17,-9))
      dev.off()
    }
    
    
  }
  
  
  # next mesoderm
  genes = c("Lefty2")
  included_colors = c("#C594BF","#1a3f52","#408DA1","#53f1fc","#DFCDE4","#B51D8D","#8DB5CE","#45d1c5")
  
  mc_f = which(mc@colors %in% included_colors)
  
  cls_f = names(mc@mc)[mc@mc %in% mc_f]
  
  for (gene in genes) {
    legc_bulk = log2(tapply(mat_n[gene,cls_f],mat@cell_metadata[cls_f,"age_group"],mean) + 1e-5)
    
    if(plot_pdf) {
      pdf(sprintf("%s/mesoderm_bulk_expr_%s.pdf",fig_dir,gene),width = 5, h = 3.5,useDingbats = F)
      plot(x = c(1:13),y = legc_bulk,xlab = "Age group",ylab = "",main = gene,cex.main = 2,pch = 19,type = "l",lwd = 4,col = "black",ylim = c(-17,-9))
      dev.off()
    } else {
      png(sprintf("%s/mesoderm_bulk_expr_%s.png",fig_dir,gene),width = 800, h = 560)
      plot(x = c(1:13),y = legc_bulk,xlab = "Age group",ylab = "",main = gene,cex.main = 2,pch = 19,type = "l",lwd = 4,col = "black",ylim = c(-17,-9))
      dev.off()
    }
    
    
  }
  
  
  
}

fig_4f = function() {
  if(!dir.exists("figs/paper_figs/fig4")) {
    dir.create("figs/paper_figs/fig4")
  }
  if(!dir.exists("figs/paper_figs/fig4/fig_4f")) {
    dir.create("figs/paper_figs/fig4/fig_4f")
  }
  plot_focals(fig_dir = "figs/paper_figs/fig4/fig_4f",focals = "Primitive streak",plot_pdf = T)
  
  
}

fig_4g = function(plot_pdf = T) {
  #show direct and ultimate fates over time and per metacell
  fig_dir = "figs/paper_figs/fig4"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  mct = scdb_mctnetwork("sing_emb_wt10")
  
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  meso_mcs = which(mc@colors %in% c("#C594BF","#DFCDE4"))
  meso_cells = names(mc@mc)[mc@mc %in% meso_mcs]
  
  mc_ag = table(mc@mc, mat@cell_metadata[names(mc@mc), "age_group"])
  #mc_ag_rn = mc_ag/rowSums(mc_ag)
  mc_ag_n = mct@mc_t_infer
  ct_ag_n = tgs_matrix_tapply(t(mc_ag_n),mc@colors,sum)
  colnames(ct_ag_n) = c(1:13)
  
  meso_fraction = colSums(ct_ag_n[c("#C594BF","#DFCDE4"),])
  
  mc_ag_rn = mct@mc_t_infer/rowSums(mct@mc_t_infer)
  rownames(mc_ag_rn) = c(1:length(mc@colors))
  colnames(mc_ag_rn)= c(1:ncol(mc_ag_rn))
  
  #mc_fate = mm_mctnetwork_get_flow_mat(mct = mct,time = "-1")
  mc_fate = mctnetwork_get_flow_mat(mct = mct,time = "-2")
  mc_fate = mc_fate[-1,-1]
  mc_fate = mc_fate/rowSums(mc_fate)
  
  meso_direct_fate = t(tgs_matrix_tapply(mc_fate[meso_mcs,],mc@colors,sum))
  rownames(meso_direct_fate) = meso_mcs
  
  cmp = cmp_commitment_vs_time(mc_list = meso_mcs,mct = mct,time_points = c(2:10),mc = mc)
  
  
  ct_nms= c("Blood progenitors","Haematoendothelial progenitors","Amnion/Chorion","ExE mesoderm","Primitive streak","Anterior Primitive Streak",
            "PGC","Cardiac mesoderm","Rostral mesoderm","Paraxial mesoderm","Lateral & intermediate mesoderm","Caudal mesoderm","Caudal epiblast")
  
  ct_rank = c(1:length(ct_nms))
  names(ct_rank) = ct_nms
  
  
  exe_ct = c("Blood progenitors","Haematoendothelial progenitors","Amnion/Chorion","ExE mesoderm")
  emb_ct = c("Primitive streak","Anterior Primitive Streak","PGC","Cardiac mesoderm","Rostral mesoderm",
             "Paraxial mesoderm","Lateral & intermediate mesoderm","Caudal mesoderm","Caudal epiblast")
  
  ct_fate = cmp$ct_fate
  ct_fate = t(t(ct_fate)*meso_fraction[as.numeric(colnames(ct_fate))])
  ct_fate = ct_fate[order(ct_rank[col_to_ct[rownames(ct_fate)]]),]
  
  exe_exc = colSums(ct_fate[col_to_ct[rownames(ct_fate)] %in% exe_ct,])
  emb_exc = colSums(ct_fate[col_to_ct[rownames(ct_fate)] %in% emb_ct,])
  max_exc = 0.3
  max_exc_neg = 0.15
  exe_exc = max_exc_neg - exe_exc
  emb_exc = max_exc - emb_exc
  
  ct_fate = rbind(exe_exc,ct_fate)
  ct_fate = rbind(ct_fate,emb_exc)
  
  bar_col = c(c("white"),rownames(ct_fate)[2:(nrow(ct_fate)-1)],c("white"))
  bar_border = c(c("white"),rep("black",nrow(ct_fate)-2),c("white"))
  if (plot_pdf) {
    pdf("figs/paper_figs/fig4/fig_4g.pdf",useDingbats = F)
  } else {
    png("figs/paper_figs/fig4/fig_4g.png")
  }
  
  barplot(ct_fate,col = bar_col,border = bar_border,axes = F)
  axis(side = 2,at = seq(0,max_exc+max_exc_neg,length.out = 4),labels = seq(-max_exc_neg,max_exc,length.out = 4))
  dev.off()
}


fig_4h = function() {
  if(!dir.exists("figs/paper_figs/fig4")) {
    dir.create("figs/paper_figs/fig4")
  }
  if(!dir.exists("figs/paper_figs/fig4/fig_4h")) {
    dir.create("figs/paper_figs/fig4/fig_4h")
  }
  plot_focals(fig_dir = "figs/paper_figs/fig4/fig_4h",focals = c("Early nascent mesoderm"),plot_pdf = T)
  
  # next follow the mc graphs for single genes
  
  mct = scdb_mctnetwork("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  confu = mctnetwork_get_flow_mat(mct, -2)
  diag(confu) = 0
  colnames(confu)[1] = "-1"

  mat = scdb_mat("sing_emb_wt10")
  mc_ag = table(mc@mc, mat@cell_metadata[names(mc@mc), "age_group"])
  mc_t = apply(mc_ag,1, function(x) sum((1:13)*x)/sum(x))
  
  
  # blood vs amnion chorion rostral mesoderm
  
  genes = c("Mesp1","Tal1","Lefty2","Hand1")
  
  incl_colors = c("#DABE99","#C594BF","#FBBE92","#DFCDE4","#cc7818","#8DB5CE","#c9a997","#8870ad","#53f1fc","#1a3f52")
  
  for(gene in genes) {
    
    plot_time_gene_mc_flow(confu = confu,mc = mc,mc_t = mc_t,gene = gene,incl_colors = incl_colors,fig_dir = "figs/paper_figs/fig4/fig_4h",max_t = 13,fig_pref = "early_nascent_mesoderm",plot_pdf = T)
  }
  
  
  
}

fig_4i = function() {
  if(!dir.exists("figs/paper_figs/fig4")) {
    dir.create("figs/paper_figs/fig4")
  }
  if(!dir.exists("figs/paper_figs/fig4/fig_4i")) {
    dir.create("figs/paper_figs/fig4/fig_4i")
  }
  plot_focals(fig_dir = "figs/paper_figs/fig4/fig_4i",focals = c("Late nascent mesoderm"),plot_pdf = T)
  
  # next follow the mc graphs for single genes
  
  mct = scdb_mctnetwork("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  confu = mctnetwork_get_flow_mat(mct, -2)
  diag(confu) = 0
  colnames(confu)[1] = "-1"
  
  mat = scdb_mat("sing_emb_wt10")
  mc_ag = table(mc@mc, mat@cell_metadata[names(mc@mc), "age_group"])
  mc_t = apply(mc_ag,1, function(x) sum((1:13)*x)/sum(x))
  
  
  # blood vs amnion chorion rostral mesoderm
  
  genes = c("Sfrp1","Dll3","Foxc1","Hand2")
  
  incl_colors = c("#DABE99","#C594BF","#FBBE92","#DFCDE4","#cc7818","#8DB5CE","#c9a997","#8870ad","#53f1fc","#1a3f52")
  
  for(gene in genes) {
    
    plot_time_gene_mc_flow(confu = confu,mc = mc,mc_t = mc_t,gene = gene,incl_colors = incl_colors,fig_dir = "figs/paper_figs/fig4/fig_4h",max_t = 13,fig_pref = "late_nascent_mesoderm",plot_pdf = T)
  }
  
  
  
}

plot_time_gene_mc_flow = function(confu, mc,mc_t, gene,incl_colors,src_y = -16, src_t=0.5,
                                  fig_dir = "figs/paper_figs", fig_pref = "epi_traj",
                                  min_t = 0.99, max_t = 10,
                                  max_y = NULL,min_y = NULL,
                                  plot_mc_ids = F,plot_pdf = F) {
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  if(plot_pdf) {
    pdf(sprintf("%s/%s.%s.pdf",fig_dir, fig_pref, gene), w=8, h=5,useDingbats = F)
  } else {
    png(sprintf("%s/%s.%s.png",fig_dir, fig_pref, gene), w=800, h=400)
  }
  
  legc = log2(mc@e_gc+1e-5)
  parent = rownames(confu)[apply(confu,2,which.max)]
  names(parent) = colnames(confu)
  
  f = (mc@colors %in% incl_colors) & (mc_t < max_t) & (mc_t > min_t)
  
  
  mc_list = as.character(which(f))
  
  parent_mc_list = parent[mc_list]
  f_arr = parent_mc_list != "-1" & parent_mc_list %in% mc_list
  
  vs = legc[gene,]
  names(vs) = as.character(1:length(vs))
  #vs["-1"] = src_y
  #mc_t["-1"] = src_t
  names(mc@colors) = as.character(1:length(mc@colors))
  if(!is.null(max_y)) {
    plot(mc_t[mc_list], legc[gene,mc_list], pch=19, col=mc@colors[mc_list],cex=1.3, xlim=c(floor(min(mc_t[parent_mc_list[f_arr]])), max(mc_t[parent_mc_list[f_arr]])),
         ylim = c(min_y,max_y),ylab=gene, xlab = "time",main = gene,cex.main = 2)
  } else {
    plot(mc_t[mc_list], legc[gene,mc_list], pch=19, col=mc@colors[mc_list],cex=1.3, xlim=c(floor(min(mc_t[parent_mc_list[f_arr]])), max(mc_t[parent_mc_list[f_arr]])), ylab=gene, xlab = "time",
         main = gene,cex.main = 2)
  }
  points(mc_t[mc_list], legc[gene,mc_list], pch=19, col=mc@colors[mc_list],cex=5)
  Arrows(mc_t[parent_mc_list[f_arr]],vs[parent_mc_list[f_arr]], mc_t[mc_list[f_arr]],vs[mc_list[f_arr]],arr.type="triangle",arr.width=0.3,arr.length=0.3, arr.adj=1,lwd = 1)
  if(plot_mc_ids) {
    text(mc_t[mc_list], legc[gene,mc_list], mc_list, cex=0.8)
  }
  
  dev.off()
}





epiblast_plot_bulk_dynamics = function(mc_list,time_points,fn,feat_genes,cluster_genes = T,n_clust = 5,highlighted_genes = NULL,
                                       seed = 12,show_lfp = F,w = 10,h = 16) {
  
  

  shades = colorRampPalette(RColorBrewer::brewer.pal(n = 9,name = "Blues")[1:9])(100)
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat_n = t(t(mat@mat[feat_genes,])/colSums(mat@mat))
  
  
  cls = names(mc@mc)[(mc@mc %in% mc_list) & (mat@cell_metadata[names(mc@mc),"age_group"] %in% time_points)]
  
  egc_t = t(tgs_matrix_tapply(mat_n[,cls],mat@cell_metadata[cls,"age_group"],mean))
  rownames(egc_t) = rownames(mat_n)
  
  legc = log2(egc_t + 1e-5)
  
  max_legc = apply(legc,1,max)
  
  legc_traj_n = (legc - log2(1e-5))/(max_legc - log2(1e-5))
  legc_traj_nn = legc_traj_n/rowSums(legc_traj_n)
  
  if(show_lfp) {
    legc_traj_n = legc[gnms,] - rowMeans(legc[gnms,])
  }
  
  if(cluster_genes) {
    km_traj = tglkmeans::TGL_kmeans(as.data.frame(legc_traj_n),k = n_clust,id_column = F,seed = seed)
    names(km_traj$cluster) = rownames(legc_traj_n)
    
    gene_mean_time = (legc_traj_nn %*% as.numeric(colnames(legc)))[,1]
    gene_mode = apply(legc_traj_nn,1,which.max)
    
    cluster_time = tapply(gene_mean_time,km_traj$cluster,mean)
    cluster_rank = rank(cluster_time)
    
    gene_rank = 1000*cluster_rank[km_traj$cluster] + 10*gene_mode + gene_mean_time
    
    
    gene_ord = order(gene_rank)
    legc_traj_n = legc_traj_n[gene_ord,]
    km_traj$cluster = km_traj$cluster[rownames(legc_traj_n)]
    
    cluster_gap_position = which(diff(km_traj$cluster) != 0)
  }
  
  
  
  breaks = seq(0,max(legc_traj_n),length.out = 101)
  
  show_rownames = F
  if(!is.null(highlighted_genes)) {
    show_rownames = T
    rownames(legc_traj_n)[!(rownames(legc_traj_n) %in% highlighted_genes)] = ""
  }
  rownames(legc_traj_n) = substr(rownames(legc_traj_n),1,10)
  
  
  if(cluster_genes) {
    pheatmap::pheatmap(legc_traj_n,cluster_rows = F,cluster_cols = F,
                       filename = fn,
                       h = h,w = w,
                       show_rownames = show_rownames,
                       show_colnames = T,
                       main = "",
                       color = shades,angle_col = 0,
                       border_color = NA,legend = F,fontsize = 20,
                       gaps_row = cluster_gap_position,breaks = breaks)
  } else {
    pheatmap::pheatmap(legc_traj_n,cluster_rows = F,cluster_cols = F,
                       filename = fn,
                       h = h,w = w,
                       show_rownames = show_rownames,
                       show_colnames = T,
                       main = "",
                       color = shades,angle_col = 0,
                       border_color = NA,legend = F,fontsize = 20,breaks = breaks)
  }
  
}





cmp_commitment_vs_time = function(mc_list,mct,time_points,mc,mc_commit = mc_list) {
  
  
  mc_ag_n_f = t(t(mct@mc_t_infer[mc_list,time_points])/colSums(mct@mc_t_infer[mc_list,time_points]))
  rownames(mc_ag_n_f) = mc_list
  colnames(mc_ag_n_f) = time_points
  
  mc_forward_m = rep(list(diag(rep(1,length(mc@colors)))),length(mct@mc_forward))
  
  for (i in 1:length(mct@mc_forward)) {
    mc_forward_m[[i]][mc_commit,] = mct@mc_forward[[i]][mc_commit,]
  }
  
  i = 13
  mc_fate_dist = diag(rep(1,length(mc@colors)))
  mc_fate_ls = list()
  mc_fate_ls[[as.character(i)]] = mc_fate_dist
  
  
  for (i in 12:1) {
    
    
    mc_fate_dist = mc_forward_m[[i]] %*% mc_fate_dist
    mc_fate_ls[[as.character(i)]] = mc_fate_dist
  }
  
  
  i = time_points[1]
  mc_fate = (mc_ag_n_f[as.character(mc_list),as.character(i)] %*% mc_fate_ls[[as.character(i)]][mc_list,])[1,]
  
  for (i in time_points[-1]) {
    fate_tmp = (mc_ag_n_f[as.character(mc_list),as.character(i)] %*% mc_fate_ls[[as.character(i)]][mc_list,])[1,]
    mc_fate = cbind(mc_fate,fate_tmp)
  }
  
  colnames(mc_fate) = time_points
  
  ct_fate = tgs_matrix_tapply(t(mc_fate),mc@colors,sum)
  ct_fate = ct_fate[rowSums(ct_fate)> 0,]
  colnames(ct_fate) = time_points
  
  return(list(mc_fate = mc_fate,ct_fate = ct_fate))
}



