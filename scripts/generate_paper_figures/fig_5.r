# fig 5 plots
#source("paper_scripts/fig3_plots.r")
source("scripts/additional_network_functions.r")
library(shape)
library(tidyverse)
library("metacell")
scdb_init("scrna_db/",force_reinit = T)

gen_fig_5_plots = function() {
  
  if(!dir.exists("figs/paper_figs/fig5")) {
    dir.create("figs/paper_figs/fig5")
  }
  if(!dir.exists("data/fig5")) {
    dir.create("data/fig5")
  }
  
  fig_5a()
  fig_5b()
  fig_5c()
  # If the two files 
  # data/fig5/best_pred_1_tf.txt
  # data/fig5/best_pred_2tfs.txt
  # are missing, please run compute_linear_model_predictions_fig_5def() before running fig_5d, fig_5e, fig_5f
  # compute_linear_model_predictions_fig_5def()
  fig_5d()
  fig_5e()
  fig_5f()
  
  # list of mesodermal transcription factors displayed in Figure 5b
  # was generated using the function
  # 
  # tf_list_mesoderm = generate_tf_mesoderm_list()
  # 
  # and saving it as tab-separated table to
  #
  # data/fig5/tf_list_mesoderm.txt

}


fig_5a = function() { 
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  ct_to_col = mc@color_key$color
  names(ct_to_col) = mc@color_key$group
  mct = scdb_mctnetwork("sing_emb_wt10")
  
  confu = mctnetwork_get_flow_mat(mct, -2)
  diag(confu) = 0
  confu = confu[-1,-1]
  ct_confu = tgs_matrix_tapply(t(confu),col_to_ct[mc@colors],sum)
  row_names = rownames(ct_confu)
  ct_confu = t(tgs_matrix_tapply(ct_confu,col_to_ct[mc@colors],sum))
  rownames(ct_confu) = row_names
  diag(ct_confu) = 0
  
  prev_ct = apply(ct_confu,2,which.max)
  prev_ct = rownames(ct_confu)[prev_ct]
  
  ct_transition = data.frame(ct1 = prev_ct,ct2 = colnames(ct_confu),stringsAsFactors = F)
  ct_transition = ct_transition[!(ct_transition$ct2 %in%c("Epiblast","Visceral endoderm","ExE endoderm","PGC")),]
  
  ct_trans_ord = c(10,22,3,12,20,11,16,18,23,8,2,15,1,5,6,19,21,17,4,13,14,7,9,24,25)
  ct_transition = ct_transition[ct_trans_ord,]
  
  
  # filter tfs that show some variance across the whole datatset
  all_tfs = read.table("data/external_data/all_tfs_list.txt",sep = "\t",stringsAsFactors = F)$x
  f = log2(apply(mc@mc_fp[all_tfs,],1,max)) > 1
  tf_f = all_tfs[f]
  
  
  diff_genes_ls = list()
  
  lfp_src_targ = c()
  
  for (i in 1:nrow(ct_transition)) {
    
    type1 = ct_to_col[ct_transition$ct1[i]]
    type2 = ct_to_col[ct_transition$ct2[i]]
    
    df_src_targ = mctnetwork_get_egc_on_cluster_transition(mct = mct,min_time = 1,max_time = 13,type1 = type1,type2 = type2)
    lfp_src_targ = cbind(lfp_src_targ,df_src_targ$lf)
    f1 = pmax(log2(df_src_targ$src + 1e-5),log2(df_src_targ$targ + 1e-5)) > -15
    f2 = abs(df_src_targ$lf) > 1
    
    diff_genes_ls[[i]] = rownames(df_src_targ)[f1 & f2]
  }
  
  diff_genes = unique(unlist(diff_genes_ls))
  feat_genes = union(diff_genes,tf_f)
  
  colnames(lfp_src_targ) = paste(ct_transition$ct1,ct_transition$ct2,sep = ",")
  rownames(lfp_src_targ) = rownames(mc@e_gc)
  lfp_src_targ = lfp_src_targ[feat_genes,]
  
  
  col_names = paste(ct_transition$ct1,ct_transition$ct2,sep = ",")
  lft = lfp_src_targ[diff_genes,]
  lft = lft[,col_names]
  
  row_max = apply(lft,1,max)
  row_max2 = apply(lft,1,FUN = function(v) {
    a = sort(-v,partial = 2)[2]
    return(-a)
  })
  
  f = (row_max > row_max2 + 0.5 & (row_max2 < 1))
  
  genes1 = rownames(lft)[f]
  genes1_ord = order(100*apply(lft[genes1,],1,which.max) + row_max[genes1])
  genes1 = genes1[genes1_ord]
  
  # plot heatmap showing genes, that are upregulated specifically in one transition between two cell types
  
  ct_col_df = mc@color_key[,c("group","color")]
  colnames(ct_col_df)= c("ct1","color1")
  ct_transition = left_join(ct_transition,ct_col_df,by = "ct1")
  colnames(ct_col_df)= c("ct2","color2")
  ct_transition = left_join(ct_transition,ct_col_df,by = "ct2")
  
  lfp_src_targ = lfp_src_targ[genes1,col_names]
  lfp_src_targ = pmax(lfp_src_targ,0)
  lfp_src_targ = pmin(lfp_src_targ,2)
  
  shades_ct_trans = colorRampPalette(RColorBrewer::brewer.pal(9,"BuPu"))(100)
  breaks_ct_trans = c(seq(0,2,length.out = 101))
  
  pdf("figs/paper_figs/fig5/fig_5a.pdf",w = 19, h = 10,useDingbats = F)
  layout(mat = matrix(c(1,2),nrow = 2,ncol = 1),heights = c(300,600))
  par(mar = c(0,3,3,3))
  par(plt = c(0,1,0,1))
  plot(x = c(1:nrow(ct_transition)),y = rep(2.5,nrow(ct_transition)),pch = 19, col = ct_transition$color1,cex = 8,ylim = c(0,4),xlim = c(0.7, nrow(ct_transition) + 0.3),axes = F)
  points(x = c(1:nrow(ct_transition)),y = rep(0.5,nrow(ct_transition)),pch = 19, col = ct_transition$color2,cex = 8)
  Arrows(x0 = c(1:nrow(ct_transition)),x1 = c(1:nrow(ct_transition)),y1 = rep(1.2,nrow(ct_transition)),y0 = rep(2,nrow(ct_transition)),lwd = 3,arr.type = "triangle",arr.width = 0.6)
  par(mar = c(3,3,0,3))
  image(t(lfp_src_targ),col = shades_ct_trans,breaks = breaks_ct_trans,axes = F)
  dev.off()
  
}

fig_5a_color_bar = function() {
  
  shades_col_bar = colorRampPalette(RColorBrewer::brewer.pal(9,"BuPu"))(101)
  vals = seq(0,2,length.out = 101)
  cols = shades_col_bar
  show_vals_ind =  c(1,51,101)
  
  
  pdf(file = "figs/paper_figs/fig5/fig_5a_color_scale.pdf",useDingbats = F)
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')
  
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  dev.off()
  
  
}


fig_5b = function() {
  
  # function creates heatmap of transcription factors in the mesoderm
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(mc@e_gc[,mc_f] + 3e-5)

  tf_f = read.table(file = "data/fig5/tf_list_mesoderm.txt",stringsAsFactors = F,sep = "\t")$x
  
  shades = colorRampPalette(RColorBrewer::brewer.pal(9,"BuPu"))(100)
  
  annotation_column = data.frame(col1 = mc@colors[mc_f])
  rownames(annotation_column) = mc_f
  col_to_col = mc@color_key$color
  names(col_to_col) = mc@color_key$color
  annotation_color = list(col1 = col_to_col)
  
  ct_n_clust = data.frame(ct = c("Early nascent mesoderm",
                                 "Late nascent mesoderm",
                                 "Caudal mesoderm",
                                 "Paraxial mesoderm",
                                 "Rostral mesoderm",
                                 "Cardiac mesoderm",
                                 "Lateral & intermediate mesoderm",
                                 "ExE mesoderm",
                                 "Allantois",
                                 "Haematoendothelial progenitors",
                                 "Amnion/Chorion"),
                          n_clust = c(5,3,2,2,4,2,3,3,1,3,1),stringsAsFactors = F)
  
  mc_clust = c()
  n_clst_tot = 0
  
  for (i in 1:nrow(ct_n_clust)) {
    ct = ct_n_clust$ct[i]
    n_clust = ct_n_clust$n_clust[i]
    mc_ct =  which(col_to_ct[mc@colors] == ct)
    
    mc_hclst = hclust(as.dist(1 - tgs_cor(legc[tf_f,as.character(mc_ct)])))
    
    mc_clust_tmp = cutree(mc_hclst,n_clust)[mc_hclst$order] + n_clst_tot
    names(mc_clust_tmp) = mc_ct[mc_hclst$order]
    mc_clust = c(mc_clust,mc_clust_tmp)
    n_clst_tot = n_clst_tot + n_clust
  }
  
  cluster_gap_position = which(diff(mc_clust) != 0)
  pheatmap::pheatmap(mat = pmin(legc[tf_f,names(mc_clust)],-10),color = shades,treeheight_row = 0,
                     cluster_cols = F,annotation_col = annotation_column,annotation_colors = annotation_color,annotation_legend = F,
                     w = 10,h = 12,filename = "figs/paper_figs/fig5/fig_5b.pdf",
                     gaps_col = cluster_gap_position,show_colnames = F,fontsize_row = 7,border_color = NA)

}



fig_5c = function() {
  
  fig_dir = "figs/paper_figs/fig5/fig_5c"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  mct = scdb_mctnetwork("sing_emb_wt10")
  mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"age_group"])
  mc_ag_n = t(t(mc_ag)/colSums(mc_ag))
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(mc@e_gc[,mc_f] + 3e-5)
  min_legc = apply(legc,1,min)
  max_legc = apply(legc,1,max)
  
  tfs = read.table(file = "data/fig5/tf_list_mesoderm.txt",sep = "\t",stringsAsFactors = F)$x
  
  for (tf in tfs) {
    
    f = (legc[tf,] - min_legc[tf])/(max_legc[tf] - min_legc[tf]) > 0.8
    
    
    mc_ff = mc_f[f]
    
    if(length(mc_ff) > 1) {
      mass_t = colSums(mc_ag_n[mc_ff,])
    } else {
      mass_t = mc_ag_n[mc_ff,]
    }
    total_mass = sum(mc_ag_n[mc_ff,])
    
    times_f = which(mass_t/total_mass > 0)
    tf_traj = rep(0,13)
    
    mc_prob_all = matrix(0,ncol = 13,nrow = ncol(mc@e_gc))
    
    for (t in times_f) {
      p_mc = rep(0,ncol(mc@e_gc))
      p_mc[mc_ff] = mc_ag_n[mc_ff,t]
      
      probs = mctnetwork_propogate_from_t(mct, t, p_mc)
      
      mc_prob_all = mc_prob_all + probs$probs
    }
    
    mc_prob_all = mc_prob_all/total_mass
    
    plot_genes_along_trajectory(mc_prob = mc_prob_all,mc = mc,fig_dir = fig_dir,genes = tf,show_sd = T,plot_pdf = T,add_gridline = T,ylim_min = -17,ylim_max = -8)
    
    
  }
  
  
  
  
}

plot_genes_along_trajectory = function(mc_prob,mc,fig_dir,genes,show_sd = T,plot_pdf = T,ylim_max = NULL,ylim_min = NULL,add_gridline = F,col = "black",tag = NULL) {
  
  
  
  mc_prob = t(t(mc_prob)/colSums(mc_prob))
  
  legc = log2(mc@e_gc + 1e-5)
  
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
  
  
  for (gene in genes) {
    df_f = df_all[df_all$gene == gene,]
    gene_name = gsub(";","_",gene)
    gene_name = gsub("/","_",gene_name)
    if(!show_sd) {
      out = ggplot(data = df_f,aes(x = age_group,y = e_gc)) +
        geom_line(aes(group = gene)) + 
        ylim(ylim_min,ylim_max) +
        ggtitle(label = gene) + 
        theme(legend.position="none") +
        theme(panel.background = element_rect(fill = "white",color = "black")) +
        theme(axis.text = element_text(size = 20),
              axis.title = element_text(size = 20),
              plot.title = element_text(size = 60,hjust = 0.5))
      
    } else {
      out = ggplot(data = df_f,aes(x = age_group,y = e_gc)) +
        xlab("Age group") +
        ylab("") +
        ylim(ylim_min,ylim_max) +
        geom_line(aes(group = gene,size = 0.5)) + 
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
      ggsave(plot = out,filename = fn,width = 10,height = 7)
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


fig_5d = function(plot_pdf = T) {
  
  # these external tables are generated by the function predict_variable_genes_lm()
  best_pred_1_tf = read.table(file = "data/fig5/best_pred_1_tf.txt",sep = "\t",stringsAsFactors = F)
  best_pred_2_tfs = read.table(file = "data/fig5/best_pred_2tfs.txt",sep = "\t",stringsAsFactors = F)
  
  #f = best_pred_2_tfs$m_rsq - best_pred_1_tf$m_rsq > 0.2 & best_pred_2_tfs$m_rsq > 0.7
  f = best_pred_2_tfs$m_rsq - best_pred_1_tf$m_rsq > 0.22 & best_pred_2_tfs$m_rsq > 0.75
  
  f2 = best_pred_1_tf$m_rsq > 0.84
  
  highlighted_tfs = best_pred_2_tfs$targ[f]
  df_plot = data.frame(targ = best_pred_2_tfs$targ[f],r_squ_1_tf = best_pred_1_tf$m_rsq[f],r_squ_2_tfs = best_pred_2_tfs$m_rsq[f],stringsAsFactors = F)
  df_plot = df_plot[order(df_plot$r_squ_2_tfs),]
  df_plot2 = data.frame(targ = best_pred_2_tfs$targ[f2],r_squ_1_tf = best_pred_1_tf$m_rsq[f2],r_squ_2_tfs = best_pred_2_tfs$m_rsq[f2],stringsAsFactors = F)
  df_plot2 = df_plot2[order(df_plot2$r_squ_1_tf),]
  
  if(plot_pdf) {
    pdf("figs/paper_figs/fig5/fig_5d.pdf",useDingbats = F)
  } else {
    png("figs/paper_figs/fig5/fig_5d.png")
  }
  
  par(mar = c(4,5,3,2))
  plot(as.numeric(best_pred_1_tf[,3]),as.numeric(best_pred_2_tfs[,4]),pch = 19,cex = 2,
       xlab = expression(paste(paste("R"^"2")," single TF")),
       ylab = expression(paste(paste("R"^"2")," two TFs")),cex.lab = 1.5)
  abline(a = 0,b = 1,lty = 2)
  
  for (i in 1:nrow(df_plot)) {
    
    segments(x0 = 0.25,x1 = df_plot$r_squ_1_tf[i],y0 = (0.95 - 0.3)*i/nrow(df_plot) + 0.3,y1 = df_plot$r_squ_2_tfs[i])
    text(x = 0.15, y = (0.95 - 0.3)*i/nrow(df_plot) + 0.3,labels = substr(df_plot$targ[i],1,10),cex = 1.5)
    
  }
  
  
  for (i in 1:nrow(df_plot2)) {
    if (i %% 2 == 0) {
      segments(x0 = 0.2 + (0.9 - 0.2)*i/nrow(df_plot2) ,x1 = df_plot2$r_squ_1_tf[i],y0 = 0.3,y1 = df_plot2$r_squ_2_tfs[i])
      text(x =0.2 + (0.9 - 0.2)*i/nrow(df_plot2), y = 0.25,labels = substr(df_plot2$targ[i],1,10),cex = 1.5)
    } else {
      segments(x0 = 0.2 + (0.9 - 0.2)*i/nrow(df_plot2) ,x1 = df_plot2$r_squ_1_tf[i],y0 = 0.25,y1 = df_plot2$r_squ_2_tfs[i])
      text(x =0.2 + (0.9 - 0.2)*i/nrow(df_plot2), y = 0.2,labels = substr(df_plot2$targ[i],1,10),cex = 1.5)
    }
    
    
  }
  
  dev.off()
  
  
  
}

fig_5e = function(plot_pdf = T) {
  
  fig_dir = "figs/paper_figs/fig5/fig_5e"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(mc@e_gc[,mc_f] + 3e-5)
  
  # these external tables are generated by the function predict_variable_genes_lm()
  best_pred_1_tf = read.table(file = "data/fig5/best_pred_1_tf.txt",sep = "\t",stringsAsFactors = F)
  best_pred_2_tfs = read.table(file = "data/fig5/best_pred_2tfs.txt",sep = "\t",stringsAsFactors = F)
  
  rownames(best_pred_1_tf) = best_pred_1_tf$targ
  rownames(best_pred_2_tfs) = best_pred_2_tfs$targ
  
  target_genes = c("Dll1","Dlk1","Bmp4","Tdgf1")
  
  for (targ in target_genes) {

    tf = best_pred_1_tf[targ,"tf1"]
    
    if(plot_pdf) {
      pdf(sprintf("%s/%s_vs_%s.pdf",fig_dir,targ,tf))
    } else {
      png(sprintf("%s/%s_vs_%s.png",fig_dir,targ,tf),w = 500,h = 500)
    }
    
    par(mar = c(6,6,2,4))
    plot(x = legc[tf,],y = legc[targ,],pch = 19,cex = 4,col = mc@colors[mc_f],xlab = tf,ylab = targ,cex.lab = 2)
    dev.off()
  }
  
}

fig_5f = function(plot_pdf = T) {
  
  fig_dir = "figs/paper_figs/fig5/fig_5f"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(mc@e_gc[,mc_f] + 3e-5)
  
  # these external tables are generated by the function predict_variable_genes_lm()
  best_pred_1_tf = read.table(file = "data/fig5/best_pred_1_tf.txt",sep = "\t",stringsAsFactors = F)
  best_pred_2_tfs = read.table(file = "data/fig5/best_pred_2tfs.txt",sep = "\t",stringsAsFactors = F)
  
  rownames(best_pred_1_tf) = best_pred_1_tf$targ
  rownames(best_pred_2_tfs) = best_pred_2_tfs$targ
  
  target_genes = c("Pcdh19","Bmp2")
  
  for (targ in target_genes) {

    tf1 = best_pred_2_tfs[targ,"tf1"]
    tf2= best_pred_2_tfs[targ,"tf2"]
    
    if(plot_pdf) {
      pdf(sprintf("%s/%s_vs_%s.pdf",fig_dir,targ,tf1))
    } else {
      png(sprintf("%s/%s_vs_%s.png",fig_dir,targ,tf1),w = 500,h = 500)
    }
    par(mar = c(6,6,2,4))
    plot(x = legc[tf1,],y = legc[targ,],pch = 19,cex = 4,col = mc@colors[mc_f],xlab = tf1,ylab = targ,cex.lab = 2)
    dev.off()
    
    if (plot_pdf) {
      pdf(sprintf("%s/%s_vs_%s.pdf",fig_dir,targ,tf2))
    } else {
      png(sprintf("%s/%s_vs_%s.png",fig_dir,targ,tf2),w = 500,h = 500)
    }
    par(mar = c(6,6,2,4))
    plot(x = legc[tf2,],y = legc[targ,],pch = 19,cex = 4,col = mc@colors[mc_f],xlab = tf2,ylab = targ,cex.lab = 2)
    dev.off()
    
    out = lm(legc[targ,] ~ legc[tf1,]+legc[tf2,])
    
    if(plot_pdf) {
      pdf(sprintf("%s/%s_vs_%s_and_%s.pdf",fig_dir,targ,tf1,tf2))
    } else {
      png(sprintf("%s/%s_vs_%s_and_%s.png",fig_dir,targ,tf1,tf2),w = 500,h = 500)
    }
    par(mar = c(6,6,2,4))
    plot(x = out$fitted.values,y = legc[targ,],pch = 19,cex = 4,col = mc@colors[mc_f],xlab = sprintf("Prediction from %s and %s",tf1,tf2),ylab = targ,cex.lab = 2)
    dev.off()
    
  }
  
  
}


compute_linear_model_predictions_fig_5def = function() {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(mc@e_gc[,mc_f] + 3e-5)
  egc = mc@e_gc[,mc_f]
  legc_n = log2(mc@e_gc[,mc_f] + 3e-5) - log2(3e-5)
  legc_all = log2(mc@e_gc + 3e-5) - log2(3e-5)
  min_legc = apply(legc,1,min)
  max_legc = apply(legc,1,max)
  diff_legc = max_legc - min_legc
  
  f1 = (max_legc - min_legc > 2) & max_legc > -12
  
  feat_genes = rownames(legc)[f1]
  bad_genes = read.table("data/external_data/sing_emb_wt10.bad_genes.txt",sep = "\t",stringsAsFactors = F)$x
  feat_genes = setdiff(feat_genes, bad_genes)
  
  tf_list = read.table("data/fig5/tf_list_mesoderm.txt",sep = "\t",stringsAsFactors = F)$x
  
  marker = setdiff(feat_genes,tf_list)
  marker = marker[diff_legc[marker] > 3]
  
  print(paste0("Number of feature genes: ", length(marker)))
  
  
  best_fit_2tf = function(targ) {
    r_sq = 0
    best_fit = c(0,0,0,0,0)
    res = data.frame()
    for(i in 1:(length(tf_list)-1)) {
      tf1 = tf_list[i]
      for(tf2 in tf_list[(i+1):length(tf_list)]) {
        m = summary(lm(legc[targ,] ~ legc[tf1,]+legc[tf2,]))
        if(m$r.squared > r_sq) {
          best_fit = c(targ,tf1, tf2, m$r.squared)
          r_sq = m$r.squared
        }
      }
    }
    return(best_fit)
  }
  
  best_pred_2_tfs = sapply(marker,FUN =  function(g) {
    
    best_fit = best_fit_2tf(g)
    
    tmp = summary(lm(legc[best_fit[1],] ~ legc[best_fit[2],]))
    r_sq_tf1 = tmp$r.squared
    tmp = summary(lm(legc[best_fit[1],] ~ legc[best_fit[3],]))
    r_sq_tf2 = tmp$r.squared
    
    best_fit = c(best_fit,c(r_sq_tf1,r_sq_tf2))
    return(best_fit)
  })
  
  best_pred_2_tfs= t(best_pred_2_tfs)
  colnames(best_pred_2_tfs) = c("targ","tf1","tf2","m_rsq","tf1_r_sq","tf2_r_sq")
  
  write.table(best_pred_2_tfs,file = "data/fig5/best_pred_2tfs.txt",sep = "\t")
  #write.table(best_pred_2_tfs_m2,file = "data/paper_data/fig5/tf_prediction/best_pred_2tfs_m2.txt",sep = "\t")
  
  best_fit_1_tf = function(targ) {
    r_sq = 0
    best_fit = c(0,0,0)
    for(i in 1:length(tf_list)) {
      tf1 = tf_list[i]
      m = summary(lm(legc[targ,] ~ legc[tf1,]))
      if(m$r.squared > r_sq) {
        best_fit = c(targ,tf1, m$r.squared)
        r_sq = m$r.squared
      }
    }
    return(best_fit)
  }
  
  best_pred_1_tf = sapply(marker,FUN =  function(g) {
    
    return(best_fit_1_tf(g))
  })
  
  best_pred_1_tf = t(best_pred_1_tf)
  colnames(best_pred_1_tf) = c("targ","tf1","m_rsq")
  
  #png("figs/paper_figs/fig5/lm_screening/r_squared_1param_vs_2_param.png",w = 1500,h = 1500)
  #plot(as.numeric(best_pred_1_tf[,3]),as.numeric(best_pred_2_tfs[,4]),pch = 19,cex = 2)
  #abline(a = 0,b = 1,lty = 2)
  #text(x = as.numeric(best_pred_1_tf[,3]),y = as.numeric(best_pred_2_tfs[,4]),labels = substr(rownames(best_pred_2_tfs),1,10),cex =1)
  #dev.off()
  
  
  write.table(best_pred_1_tf,file = "data/fig5/best_pred_1_tf.txt",sep = "\t")
  
  
  
  
}


generate_tf_mesoderm_list = function() {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(mc@e_gc[,mc_f] + 3e-5)
  min_legc = apply(legc,1,min)
  max_legc = apply(legc,1,max)
  
  f1 = (max_legc - min_legc > 2) & max_legc > -12
  
  feat_genes = rownames(legc)[f1]
  bad_genes = read.table("data/external_data/sing_emb_wt10.bad_genes.txt",sep = "\t",stringsAsFactors = F)$x
  feat_genes = setdiff(feat_genes, bad_genes)
  ignore_tfs = c("Ccnd1","Naca","Hmga2","Psip1","Gtf2f2","Sall4","Basp1","Smad1","Actn1","Mycn")
  
  tf_f = read.table("data/external_data/all_tfs_list.txt",sep = "\t",stringsAsFactors = F)
  tf_f = tf_f$x
  
  tf_f = intersect(tf_f,feat_genes)
  tf_f = setdiff(tf_f,ignore_tfs)
  
  if (0) {
    write.table(x = tf_f,file = "data/fig5/tf_list_mesoderm.txt",sep = "\t")
  }

  
  return(tf_f)
}