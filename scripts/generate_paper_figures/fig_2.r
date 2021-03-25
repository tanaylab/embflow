source("scripts/generate_paper_figures/plot_network.r")

gen_fig_2_plots = function() {
  tgconfig::override_params("config/sing_emb.yaml","metacell")
  
  if(!dir.exists("figs/paper_figs/fig2")) {
    dir.create("figs/paper_figs/fig2")
  }
  
  fig_2a()
  fig_2b()
  fig_2c()
  fig_2d()
  fig_2d_legend()
  fig_2e()
  
}



fig_2a = function(plot_pdf = F) {
  
  
  edge_w = 1
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mc2d = scdb_mc2d("sing_emb_wt10_recolored")
  
  if(plot_pdf) {
    
    fn= "figs/paper_figs/fig2/fig_2a.pdf"
    
    res = 72
    mcp_2d_height = 1000
    mcp_2d_width = 1000
    mcp_2d_cex = 0.8
    
    pdf(file = fn,width = mcp_2d_width/res,height = mcp_2d_height/res,useDingbats = F)
    cols = mc@colors
    cols[is.na(cols)] = "gray"
    plot(mc2d@sc_x, mc2d@sc_y, pch=19, col=cols[mc@mc[names(mc2d@sc_x)]],axes = F,xlab = "",ylab = "")
    fr = mc2d@graph$mc1
    to = mc2d@graph$mc2
    
    dx = mc2d@mc_x[fr]-mc2d@mc_x[to]
    dy = mc2d@mc_y[fr]-mc2d@mc_y[to]
    segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to], 
             lwd=edge_w)
    
    points(mc2d@mc_x, mc2d@mc_y, cex= 3*mcp_2d_cex, col="black", pch=21, bg=cols)
    text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=mcp_2d_cex)
    
    dev.off()
  } else {
    
    fn= "figs/paper_figs/fig2/fig_2a.png"
    
    scl = 3
    mcp_2d_height = 1000
    mcp_2d_width = 1000
    mcp_2d_cex = 0.8
    
    png(file = fn,width = mcp_2d_width*scl,height = mcp_2d_height*scl)
    cols = mc@colors
    cols[is.na(cols)] = "gray"
    plot(mc2d@sc_x, mc2d@sc_y, pch=19, col=cols[mc@mc[names(mc2d@sc_x)]],axes = F,xlab = "",ylab = "",cex = scl)
    fr = mc2d@graph$mc1
    to = mc2d@graph$mc2
    
    dx = mc2d@mc_x[fr]-mc2d@mc_x[to]
    dy = mc2d@mc_y[fr]-mc2d@mc_y[to]
    segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to], 
             lwd=edge_w*scl)
    
    points(mc2d@mc_x, mc2d@mc_y, cex= 3*mcp_2d_cex*scl, col="black", pch=21, bg=cols,lwd = scl)
    text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=mcp_2d_cex*scl)
    
    dev.off()
    
    
  }
  
  
  
  
}


fig_2b = function(plot_pdf = T) {
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  emb_counts = read.table(file = "data/external_data/counted_embryos/wt10.cell_counts.txt",sep = "\t",stringsAsFactors = F)
  emb_time = unique(mat@cell_metadata[names(mc@mc),c("embryo","developmental_time")])
  
  emb_to_time = emb_time$developmental_time
  names(emb_to_time) = emb_time$embryo
  fit_count = lm(log2(emb_counts$cell_count) ~ emb_to_time[emb_counts$embryo])
  
  y_fit = 2^fit_count$fitted.values
  
  if(plot_pdf) {
    pdf(file = "figs/paper_figs/fig2/fig_2b.pdf",useDingbats = F)
  } else {
    png(file = "figs/paper_figs/fig2/fig_2b.png")
  }
  plot(emb_to_time[emb_counts$embryo],emb_counts$cell_count,pch = 19,log = "y",cex = 2,xlab = "Time",ylab = "Number of cells")
  lines(x = emb_to_time[emb_counts$embryo],y = y_fit)
  dev.off()
  
}


fig_2c = function(plot_pdf = FALSE) {
  
  mat_id = "sing_emb_wt10"
  mc_id = "sing_emb_wt10_recolored"
  dir_name = "figs/paper_figs/fig2"
  m_0 = 0.01
  s_0 = 0.005
  
  m = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  mc2d = scdb_mc2d("sing_emb_wt10_recolored")
  
  col_to_rank = c(1:nrow(mc@color_key))
  names(col_to_rank) = mc@color_key$color
  
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  m_genes = c("Mki67","Cenpf","Top2a","Smc4;SMC4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr;IHABP","Cenpa;Cenp-a","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgol2","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")
  s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")	# Npm1
  
  mc_mean_age = tapply(m@cell_metadata[names(mc@mc),"age_group"],mc@mc,mean)
  
  s_genes = intersect(rownames(mc@mc_fp), s_genes)
  m_genes = intersect(rownames(mc@mc_fp), m_genes)
  
  tot  = colSums(m@mat)
  s_tot = colSums(m@mat[s_genes,])
  m_tot = colSums(m@mat[m_genes,])
  
  s_score = s_tot/tot
  m_score = m_tot/tot
  
  f = (m_score < m_0 * (1- s_score/s_0))
  
  mc_cc_tab = table(mc@mc, f[names(mc@mc)])
  mc_cc = 1+floor(99*mc_cc_tab[,2]/rowSums(mc_cc_tab))
  
  
  # 2d projection and barplot
  shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)
  
  if(plot_pdf) {
    pdf(sprintf("%s/fig_2c_2d_projection.pdf", dir_name), w=16, h=16)
    plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=0.8, col=ifelse(f[names(mc2d@sc_x)], "black", "lightgray"),axes = F,xlab = "",ylab = "")
    points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=3, bg=shades[mc_cc])
    dev.off()
  } else {
    png(sprintf("%s/fig_2c_2d_projection.png", dir_name), w=1600, h=1600)
    plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=1, col=ifelse(f[names(mc2d@sc_x)], "black", "lightgray"),axes= F,xlab = "",ylab = "")
    points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=5, bg=shades[mc_cc])
    dev.off()
  }
  
  # next the color legend
  
  shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(101)
  
  vals = seq(0,1,length.out = 101)
  cols = shades
  show_vals_ind =  c(1,51,101)
  
  
  pdf(file = "figs/paper_figs/fig2/fig_2c_color_legend.pdf",useDingbats = F)
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')
  
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  dev.off()
  
  # boxplot per cell type of the cell cycle score
  
  cc_score_all = m_score + s_score
  
  all_cls = intersect(colnames(m@mat),names(mc@mc))

  mc_col_m = rep("gray",length(mc@colors))
  f = mc@colors %in% c("#F6BFCB","#7F6874","#0F4A9C")
  mc_col_m[f] = mc@colors[f]
  
  col_order = c("gray","#F6BFCB","#7F6874","#0F4A9C")
  
  
  
  cc_score = split(cc_score_all[all_cls],f = mc_col_m[mc@mc[all_cls]])
  cc_score = cc_score[col_order]
  
  pdf(sprintf("%s/fig_2c_boxplot.pdf",dir_name),useDingbats = F)
  boxplot(cc_score,pch = 19,cex = 0.5,col = names(cc_score),xaxt = 'n')
  dev.off()
  
}

fig_2d = function() {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  network_color_ord = mc@color_key$color
  
  net_id = "sing_emb_wt10"
  
  # this function is sourced from scripts/plot_network.r
  mm_mctnetwork_plot_net(mct_id = net_id,fn = "figs/paper_figs/fig2/fig_2d_flow_chart.png",w = 3500, h = 4900,
                         dx_back = 0,colors_ordered=network_color_ord,plot_pdf = F,show_axes = F,show_over_under_flow = F,mc_cex = 1,max_lwd = 15,edge_w_scale = 2e-4,
                         plot_mc_ids = F)
}

fig_2d_legend = function() {
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  
  pdf("figs/paper_figs/fig2/fig_2d_cell_type_legend.pdf",h = 8,w = 4,useDingbats = F)
  plot.new()
  legend(x = "topleft",legend = mc@color_key$group,pch = 19, col = mc@color_key$color,pt.cex = 2,cex = 1,bty = 'n')
  dev.off()
  
  
}

fig_2e = function() {
  
  fig_dir = "figs/paper_figs/fig2"
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  
  shades = colorRampPalette(RColorBrewer::brewer.pal(9,"BuPu"))(1000)
  
  network_color_ord = mc@color_key$color
  
  mc_rank = mctnetwork_mc_rank_from_color_ord(mct_id = "sing_emb_wt10",colors_ordered =  network_color_ord)
  mc_ord = order(mc_rank)

  lfp = log2(mc@e_gc + 1e-5)
  
  marker_genes = c("Utf1","Pou3f1","Sox3","Irx3","Tcfap2a","Six3","Dppa3","Grsf1","Aldh1a2","Tbx6","Eomes","Snai1","Mesp1","Tcf15","Twist1","Smarcd3","Nkx2-5","Tagln","Hand1","Tbx4","Etv2",
                   "Lmo2","Tal1","Gata1","Noto","Foxa2","Cer1","Sox17","Ttr")
  
  lfp = lfp- rowMeans(lfp)
  lfp = lfp[marker_genes,]
  lfp = lfp[,rev(mc_ord)]
  
  lfp = pmax(lfp,0)
  lfp = pmin(lfp,6)
  
  lfp_n = lfp/rowSums(lfp)
  
  cluster_gap_position = c(2,3,4,5,6,7,10,13,15,17,20,24,25,28)
  
  annotation_row = data.frame(ct = mc@colors)
  rownames(annotation_row) = c(1:length(mc@colors))
  col_to_col = unique(mc@colors)
  names(col_to_col) = col_to_col
  annotation_color = list(ct = col_to_col)
  
  fn = sprintf("%s/fig_2e.pdf",fig_dir)
  
  pheatmap::pheatmap(mat = t(lfp),cluster_cols = F,cluster_rows = F,color = shades,filename = fn,show_rownames = F,
                     gaps_col = cluster_gap_position,
                     column_gap = unit(100,"mm"),
                     annotation_row = annotation_row,annotation_colors = annotation_color,
                     annotation_legend = F,annotation_names_row = F,w = 6.4, h = 10.4)
    
  
  
    
}

