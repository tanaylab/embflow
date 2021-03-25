library("metacell")
scdb_init("scrna_db/",force_reinit = T)
library(gridExtra)



gen_fig_1_plots = function() {
  
  if(!dir.exists("figs/paper_figs")) {
    dir.create("figs/paper_figs")
  }
  dir_name = "figs/paper_figs/fig1"
  
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  fig1_b()
  fig1_cde()
  fig1_f()
  fig1_g_mc_time_distributions()
  fig1_g_heatmap(plot_pdf = T)
  fig1_h()
  
}



fig1_cde = function() {
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  emb_ranks = unique(mat@cell_metadata[names(mc@mc),c("embryo","transcriptional_rank","ref_gastru_atlas_rank","ref_gastru_atlas_age","morphology_rank","area","sex")])
  #emb_ranks$color = as.character(ifelse(emb_ranks$sex == "female","coral3","cornflowerblue"))
  emb_ranks$color = as.character(ifelse(emb_ranks$sex == "female","#CB181D","cornflowerblue"))

  
  dir_name = "figs/paper_figs/fig1"
  
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  cex.lab = 2
  cex.axis = 2
  cex.main = 2
  margins = c(5,5,5,5)
  cex = 3
  
  
  xlabel = ""
  ylabel = ""
  main_label = ""
  #xlabel = "intrinsic rank"
  #ylabel = "atlas rank"
  #main_label = "Atlas vs Intrinsic Rank"
  
  pdf(sprintf("%s/fig1_d.pdf",dir_name),useDingbats = F)
  #png(sprintf("%s/atlas_rank_vs_intrinsic_rank.png",dir_name))
  par(mar = margins)
  plot(x = emb_ranks$transcriptional_rank, y = emb_ranks$ref_gastru_atlas_rank,pch = 19,cex.lab = cex.lab,cex.axis = cex.axis,
       xlab = xlabel,ylab = ylabel,col = emb_ranks$color,
       main = main_label,cex.main = cex.main,cex = cex)
  dev.off()
  
  
  xlabel = ""
  ylabel = ""
  main_label = ""  
  #xlabel = "intrinsic rank"
  #ylabel = "atlas mean age"
  #main_label = "Atlas Age vs Intrinsic Rank"
  
  if(0) {
    pdf(sprintf("%s/atlas_age_vs_intrinsic_rank.pdf",dir_name),useDingbats = F)
    #png(sprintf("%s/atlas_age_vs_intrinsic_rank.png",dir_name))
    par(mar = margins)
    plot(x = emb_ranks$transcriptional_rank, y = emb_ranks$ref_gastru_atlas_age,pch = 19,cex.lab = cex.lab,cex.axis = cex.axis,
         xlab = xlabel,ylab = ylabel,col = emb_ranks$color,
         main = main_label,cex.main = cex.main,cex = cex)
    dev.off()
  }

  f = !is.na(emb_ranks$morphology_rank)
  emb_ranks = emb_ranks[f,]
  emb_ranks$morphology_rank = rank(emb_ranks$morphology_rank)
  
  
  xlabel = ""
  ylabel = ""
  main_label = ""  
  #xlabel = "intrinsic rank"
  #ylabel = "morphological rank"
  #main_label = "Morphological vs Intrinsic Rank"
  
  pdf(sprintf("%s/fig1_c.pdf",dir_name),useDingbats = F)
  #png(sprintf("%s/morphology_rank_vs_intrinsic_rank.png",dir_name))
  par(mar = margins)
  plot(x = emb_ranks$transcriptional_rank, y = emb_ranks$morphology_rank,pch = 19,cex.lab = cex.lab,cex.axis = cex.axis,
       xlab = xlabel,ylab = ylabel,col = emb_ranks$color,
       main = main_label,cex.main = cex.main,cex = cex)
  dev.off()
  
  f = !is.na(emb_ranks$area)
  emb_ranks = emb_ranks[f,]

  if(0) {
    xlabel = ""
    ylabel = ""
    main_label = ""  
    #xlabel = "size rank"
    #ylabel = "morphological rank"
    #main_label = "Morphology vs Size"
    
    pdf(sprintf("%s/morphology_rank_size_rank.pdf",dir_name),useDingbats = F)
    #png(sprintf("%s/morphology_rank_size_rank.png",dir_name))
    par(mar = margins)
    plot(x = rank(emb_ranks$area), y = emb_ranks$morphology_rank,pch = 19,cex.lab = cex.lab,cex.axis = cex.axis,
         xlab = xlabel,ylab = ylabel,col = emb_ranks$color,
         main = main_label,cex.main = cex.main,cex = cex)
    dev.off()
    
    xlabel = ""
    ylabel = ""
    main_label = ""  
    #xlabel = "atlas mean age"
    #ylabel = "log2(area)"
    #main_label = "Atlas Age vs Size"
    
    
    pdf(sprintf("%s/atlas_age_vs_size.pdf",dir_name),useDingbats = F)
    #png(sprintf("%s/atlas_age_vs_size.png",dir_name))
    par(mar = margins)
    plot(x = emb_ranks$ref_gastru_atlas_age, y = log2(emb_ranks$area),pch = 19,cex.lab = cex.lab,cex.axis = cex.axis,
         xlab = xlabel,ylab = ylabel,col = emb_ranks$color,
         main = main_label,cex.main = cex.main,cex = cex)
    dev.off()
    
    

  }
 
  xlabel = ""
  ylabel = ""
  main_label = ""  
  #xlabel = "intrinsic rank"
  #ylabel = "log2(area)"
  #main_label = "Transcriptional Rank vs Size"
  
  pdf(sprintf("%s/fig1_e.pdf",dir_name),useDingbats = F)
  #png(sprintf("%s/intrinsic_rank_vs_size.png",dir_name))
  par(mar = margins)
  plot(x = emb_ranks$transcriptional_rank, y = log2(emb_ranks$area),pch = 19,cex.lab = cex.lab,cex.axis = cex.axis,
       xlab = xlabel,ylab = ylabel,col = emb_ranks$color,
       main = main_label,cex.main = cex.main,cex = cex)
  dev.off()
  
  
}







fig1_b = function() {
  
  dir_name = "figs/paper_figs/fig1"
  
  shades = rev(colorRampPalette(RColorBrewer::brewer.pal(9,"RdGy"),bias = 0.5)(1000))

  
  intrinsic_ranks = read.table("data/intrinsic_temporal_ranking/intrinsic_ranking_embryo_final_order_second_iteration.txt",sep= "\t",stringsAsFactors = F)
  emb_mat = read.table("data//intrinsic_temporal_ranking/intrinsic_ranking_similarity_mat_second_iteration.txt",sep= "\t",stringsAsFactors = F)
  colnames(emb_mat) = rownames(emb_mat)
  
  emb_mat = emb_mat[intrinsic_ranks$embryo,intrinsic_ranks$embryo]
  emb_mat = as.matrix(emb_mat)
  emb_mat = pmin(emb_mat,0.05)
  
  colnames(emb_mat) = rep(" ",ncol(emb_mat))
  rownames(emb_mat) = rep(" ",ncol(emb_mat))
  colnames(emb_mat)[c(20,40,60,80,100,120,140)] = c(20,40,60,80,100,120,140)
  rownames(emb_mat)[c(20,40,60,80,100,120,140)] = c(20,40,60,80,100,120,140)
  
  pheatmap::pheatmap(emb_mat[c(nrow(emb_mat):1),],cluster_rows = F,cluster_cols = F,
                     show_rownames = T,show_colnames = F,
                     col = shades,w = 5,h = 4.5,border_color = NA,
                     filename = sprintf("%s/fig1_b.pdf",dir_name))
  
}


fig1_g_mc_time_distributions = function() {
  
  # six metacells
  # epiblast
  # nascent mesoderm
  # definitive endoderm
  # neural ectoderm
  # caudal mesoderm
  # cranial mesoderm
  
  fig_dir = "figs/paper_figs/fig1/fig1_g"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  
  # color used in barplot
  gene_col = "#0570B0"
  
  #mcs = c(272,331,449,49,11,397)
  mcs = c(272,449,49,397)
  
  mc_col = mc@colors[mcs]
  
  col_palette = rev(RColorBrewer::brewer.pal(9,"PuBu"))
  col_palette[9] = "white"
  shades = rev(colorRampPalette(col_palette)(1000))
  #shades = shades[10:1000]

  
  mc_emb = table(mc@mc,mat@cell_metadata[names(mc@mc),"transcriptional_rank"])
  mc_emb_c = t(t(mc_emb)/colSums(mc_emb))
  mc_emb_cn = mc_emb_c/rowSums(mc_emb_c)  
  
  mc_mean_rank = (mc_emb_cn %*% as.numeric(colnames(mc_emb_cn)))[,1]
  
  mc_emb_cn = mc_emb_cn[order(mc_mean_rank),]
  
  mc_emb_cn = mc_emb_cn[as.character(mcs),]
  
  dens_est = apply(mc_emb_cn,1,FUN = function(mc_weight) {
    a = density(x = as.numeric(colnames(mc_emb_cn)),weights = mc_weight,from = min(as.numeric(colnames(mc_emb))),
                to = max(as.numeric(colnames(mc_emb))), n =500)
    
    a$y = a$y/sum(a$y)
    return(a)
  })
  
  legc = log2(mc@e_gc+1e-5)
  lfp = log2(mc@mc_fp)
  
  for (i in 1:nrow(mc_emb_cn)) {
    
    m = rownames(mc_emb_cn)[i]
    
    dens = dens_est[[i]]
    
    
    genes = tail(names(sort(log2(mc@mc_fp[,m]))),3)
    #genes[genes == "AK145379;H19"] = "H19"
    
    gene_names = genes
    gene_names[gene_names == "AK145379;H19"] = "H19"
    names(gene_names) = genes
    
    lfp_f = lfp[genes, m]
    names(lfp_f) = gene_names
    
    pdf(file = sprintf("%s/transcriptional_rank_distribution_mc_%s.pdf",fig_dir,m), w = 28, h = 7,useDingbats = F)
    layout(mat = matrix(c(1,2),nrow = 1,ncol = 2),widths = c(8,10))
    par(mar = c(2,20,2,1.5))
    barplot(lfp_f, las=2,col = gene_col,horiz = T,cex.names = 6,xlab = "lfp")
    par(mar = c(2,1,2,1.5))
    plot(x = dens$x,y = dens$y,type = "l",xlab = "Embryo transcriptional rank",ylab = "n",
         lwd = 10,col = mc@colors[as.numeric(m)])
    dev.off()

  }
  
}


fig1_g_heatmap = function(plot_pdf = F) {
  
  mcs = c(272,449,49,397)
  
  embryos = c("EXE3_m1_e3","190313_W_wt05","190710_C3","190223_R_b2")
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  
  emb_vs_rank = unique(mat@cell_metadata[names(mc@mc),c("embryo","transcriptional_rank")])
  emb_ranks = emb_vs_rank$transcriptional_rank[emb_vs_rank$embryo %in% embryos]
  
  col_palette = rev(RColorBrewer::brewer.pal(9,"PuBu"))
  col_palette[9] = "white"
  shades = rev(colorRampPalette(col_palette,bias = 0.7)(1000))
  #shades = shades[10:1000]
  
  mc_emb = table(mc@mc,mat@cell_metadata[names(mc@mc),"transcriptional_rank"])
  mc_emb_c = t(t(mc_emb)/colSums(mc_emb))
  mc_emb_cn = mc_emb_c/rowSums(mc_emb_c)  
  
  mc_mean_rank = (mc_emb_cn %*% as.numeric(colnames(mc_emb_cn)))[,1]
  mc_rank = rank(mc_mean_rank)
  
  mc_emb_cn = mc_emb_cn[order(mc_mean_rank),]
  
  dens_est = apply(mc_emb_cn,1,FUN = function(mc_weight) {
    a = density(x = as.numeric(colnames(mc_emb_cn)),weights = mc_weight,from = min(as.numeric(colnames(mc_emb))),
                to = max(as.numeric(colnames(mc_emb))), n =500)
    return(a$y)
  })

  if (plot_pdf) {
    pdf(file = "figs/paper_figs/fig1/fig1_g/heatmap_mc_vs_embryo_ranks.pdf", w = 10,h = 11,useDingbats = F)
    par(mar = c(5.1, 8, 4.1, 2.1))
    image(x = dens_est,col = shades,axes = F)
    axis(side = 1,at = c(emb_ranks)/153,labels = emb_ranks,lwd = 0,lwd.ticks = 1,cex.axis = 1.7)
    #axis(side = 2,labels = rownames(mc_emb_cn),tick = F,at = seq(0,1,length.out = nrow(mc_emb_cn)))
    #mtext(text = rownames(mc_emb_cn),at = seq(0,1,length.out = 6),side = 2,las = 2)
    abline(h = 1 + seq(0,1,length.out = nrow(mc_emb_cn))[2]/2)
    abline(h = 0 - seq(0,1,length.out = nrow(mc_emb_cn))[2]/2)
    abline(v = 1)
    abline(v = 0)
    axis(side = 2, at = mc_rank[mcs]/length(mc_rank),labels = paste0(rep("MC ",length(mcs)),mcs),lwd = 0,lwd.ticks = 1,las = 2,cex.axis = 1.7)
    
    for(e_r in emb_ranks) {
      abline(v = e_r/max(emb_vs_rank$transcriptional_rank),col = scales::alpha("grey30",0.3),lty = 2,lwd = 2)
    }
    
    dev.off()
  } else {
    fig_scale = 4
    png(filename = "figs/paper_figs/fig1/fig1_g/heatmap_mc_vs_embryo_ranks.png", w = 500*fig_scale,h = 550*fig_scale)
    par(mar = fig_scale*c(5.1, 6, 4.1, 2.1))
    image(x = dens_est,col = shades,axes = F)
    #axis(side = 1,at = c(emb_ranks)/153,labels = emb_ranks,lwd = 0,lwd.ticks = 1,cex.axis = fig_scale*1.7)
    #axis(side = 2,labels = rownames(mc_emb_cn),tick = F,at = seq(0,1,length.out = nrow(mc_emb_cn)))
    #mtext(text = rownames(mc_emb_cn),at = seq(0,1,length.out = 6),side = 2,las = 2)
    abline(h = 1 + seq(0,1,length.out = nrow(mc_emb_cn))[2]/2)
    abline(h = 0 - seq(0,1,length.out = nrow(mc_emb_cn))[2]/2)
    abline(v = 1)
    abline(v = 0)
    axis(side = 2, at = mc_rank[mcs]/length(mc_rank),labels = paste0(rep("MC ",length(mcs)),mcs),lwd = 0,lwd.ticks = fig_scale*1,las = 2,cex.axis = fig_scale*1.7)
    
    for(e_r in emb_ranks) {
      abline(v = e_r/max(emb_vs_rank$transcriptional_rank),col = scales::alpha("grey30",0.3),lty = 2,lwd = fig_scale*2)
    }
    
    dev.off()
  }
  
  # next plot color scale
  fig_dir = "figs/paper_figs/fig1/fig1_g"
  max_dens = max(dens_est)
  min_dens = 0
  
  
  # plot color scale
  shades = rev(colorRampPalette(col_palette,bias = 0.7)(1001))
  cols = shades
  vals = seq(min_dens,max_dens,length.out = 1001)
  
  show_val_ind = floor(1 +  c(0,0.01,0.02,0.026)*1000/max_dens)
  vals[show_val_ind] = c(0,0.01,0.02,0.026)
  
  pdf(file = "figs/paper_figs/fig1/fig1_g/heatmap_color_scale.pdf",useDingbats = F)
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')
  
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  dev.off()
  
  
  
}




fig1_h = function() {
  
  if(!dir.exists("figs/paper_figs/fig1/fig1_h")) {
    dir.create("figs/paper_figs/fig1/fig1_h")
  }
  embryos = c("EXE3_m1_e3","190313_W_wt05","190710_C3","190223_R_b2")
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  emb_cells = table(mat@cell_metadata[names(mc@mc),"embryo"])
  
  emb_vs_age = unique(mat@cell_metadata[names(mc@mc),c("embryo","transcriptional_rank","age_group")])
  emb_vs_age$ncells = emb_cells[emb_vs_age$embryo]
  
  emb_vs_age = emb_vs_age[emb_vs_age$embryo %in% embryos,]
  emb_vs_age = emb_vs_age[order(emb_vs_age$transcriptional_rank),]
  
  emb_ranks = emb_vs_age$transcriptional_rank
  
  col_palette = rev(RColorBrewer::brewer.pal(9,"PuBu"))
  col_palette[9] = "white"
  shades = rev(colorRampPalette(col_palette)(1000))
  #shades = shades[10:1000]
  
  mc_emb = table(mc@mc,mat@cell_metadata[names(mc@mc),"transcriptional_rank"])
  mc_emb = mc_emb[,emb_vs_age$transcriptional_rank]
  mc_emb_n = t(t(mc_emb)/colSums(mc_emb))
  
  mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"age_group"])
  
  mc_ag_c = t(t(mc_ag)/colSums(mc_ag))
  mc_ag_cn = mc_ag_c/rowSums(mc_ag_c)
  
  mc_mean_age = (mc_ag_cn %*% c(1:13))[,1]
  
  mc_mean_age = tapply(mat@cell_metadata[names(mc@mc),"age_group"],mc@mc,mean)
  
  mc_gotg_age = tapply(mat@cell_metadata[names(mc@mc),"ref_gastru_atlas_age"],mc@mc,mean)
  
  dens_est = apply(mc_emb_n,2,FUN = function(mc_weight) {
    a = density(x = mc_gotg_age,weights = mc_weight,from = 6.7,to = 8.2, n =500,bw = 0.05)
    a$y = a$y/sum(a$y)
    return(a)
  })
  

  for(i in 1:ncol(mc_emb_n)) {
    emb = colnames(mc_emb_n)[i]
    
    dens = dens_est[[i]]
    
    pdf(file = sprintf("figs/paper_figs/fig1/fig1_h/emb_%s.pdf",emb),useDingbats = F)
    plot(x = dens$x,y = dens$y,type = "l",xlab = "",ylab = "",
         lwd = 10)
    dev.off()
    
  }
 
   
}


fig1_f = function() {
  
  df_germ_layer = data.frame(germ_layer = c("endo","meso","ecto"),shade = c("orange3","#0570B0","#CB181D"),stringsAsFactors = F)
  
  for (i in 1:nrow(df_germ_layer)) {
    mat_name = sprintf("sing_emb_wt10_%s",df_germ_layer[i,"germ_layer"])
    mc_id = mat_name
    mc2d_id = mat_name
    mat_id = mat_name
    
    sc_color = df_germ_layer[i,"shade"]
    bg_col="grey70"
    
    dir_name = "figs/paper_figs/fig1/fig1_f"
    if(!dir.exists(dir_name)) {
      dir.create(dir_name)
    }
    mc = scdb_mc(mc_id)
    mc2d = scdb_mc2d(mc2d_id)
    mat = scdb_mat(mat_id)
    
    png(sprintf("%s/%s.2d_all.png",dir_name,mc2d_id),w = 250,h = 2600)
    layout(matrix(1:13,ncol = 1,nrow = 13),heights = rep(200,13))
    par(bty = "n")
    par(mar = c(0.4,0.4,0.4,0.4))
    for (i in 1:13) {
      f = mat@cell_metadata[names(mc@mc),"age_group"] == i
      
      plot(mc2d@sc_x[names(mc@mc)], mc2d@sc_y[names(mc@mc)], pch=19, col=bg_col,axes = F,xlab = "",ylab = "",cex = 0.9)
      points(mc2d@sc_x[names(mc@mc)[f]],mc2d@sc_y[names(mc@mc)[f]],pch = 19,col = sc_color,cex = 0.8)
      
      
    }
    dev.off()

    for (i in 1:13) {
      f = mat@cell_metadata[names(mc@mc),"age_group"] == i
      
      png(sprintf("%s/%s.2d_age%d.png",dir_name,mc2d_id,i))
      par(bty = "n")
      plot(mc2d@sc_x[names(mc@mc)], mc2d@sc_y[names(mc@mc)], pch=19, col=bg_col,axes = F,xlab = "",ylab = "",cex = 0.9)
      points(mc2d@sc_x[names(mc@mc)[f]],mc2d@sc_y[names(mc@mc)[f]],pch = 19,col = sc_color,cex = 0.8)
      dev.off()
      
    }
    
  }
  

}
