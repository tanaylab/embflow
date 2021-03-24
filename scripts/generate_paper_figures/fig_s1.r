# s1 panels
#source("paper_scripts/calculate_embryo_time.r")
library("Matrix")

gen_fig_s1_plots = function() {
  
  if(!dir.exists("figs/paper_figs/fig_s1")) {
    dir.create("figs/paper_figs/fig_s1")
  }
  
  fig_s1_a_number_of_umi_per_cell(T)
  fig_s1_a_number_of_cells_per_embryo()
  fig_s1_a_number_of_cls_per_age_group()
  
  fig_s1c()
  fig_s1d(T)
  fig_s1e()
  fig_s1f(T)
  
  fig_s1g()
  fig_s1h()
  
  fig_s1i()
  
}


fig_s1d = function(plot_pdf = T) {
  
  fig_dir = "figs/paper_figs/fig_s1"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mat_id = "sing_emb_wt10"
  mc_id = "sing_emb_wt10_recolored"

  mat = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  gene_intervals = read.table(file = "data/external_data/gene_intervals_mm9.txt",sep = "\t",h = T,stringsAsFactors = F)
  chrom_gene_table = unique(gene_intervals[,c("chrom","gene_name")])
  genes_y = chrom_gene_table$gene_name[chrom_gene_table$chrom == "chrY"]
  
  mat_n = t(t(mat@mat)/colSums(mat@mat))
  mat_ig_n = t(t(mat@ignore_gmat)/colSums(mat@mat[,colnames(mat@ignore_gmat)]))
  
  y_genes_per_embryo = colSums(mat_n[genes_y,names(mc@mc)])
  y_genes_per_embryo = tapply(X = y_genes_per_embryo,
                              INDEX = mat@cell_metadata[names(mc@mc),"embryo"],
                              FUN = mean)
  
  Xist_per_embryo = tapply(X = mat_ig_n["Xist",names(mc@mc)],
                           INDEX = mat@cell_metadata[names(mc@mc),"embryo"],
                           FUN = mean)
  
  y_umis_per_cell = colSums(mat@mat[genes_y,names(mc@mc)])
  xist_per_cell = mat@ignore_gmat["Xist",names(mc@mc)]
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }

  embryo_sex = rep("female",length(Xist_per_embryo))
  embryo_sex[y_genes_per_embryo > 5e-5] = "male"
  df_emb_sex = data.frame(embryo = names(y_genes_per_embryo), sex = embryo_sex)
  
  df_emb_sex$color = as.character(ifelse(df_emb_sex$sex == "female","#CB181D","cornflowerblue"))
  
  # Next plot per embryo
  if(plot_pdf) {
    pdf(sprintf("%s/fig_s1d.pdf",fig_dir),useDingbats = F)
  } else {
    png(sprintf("%s/fig_s1d.png",fig_dir))
  }
  plot(x = Xist_per_embryo,y = y_genes_per_embryo[names(Xist_per_embryo)],pch = 19,
       xlab = "Xist UMIs per embryo", ylab = "Y chromosome UMIs per embryo",col = df_emb_sex$color,cex = 2,cex.lab = 1.5)
  dev.off()
  
  
}


fig_s1_a_number_of_umi_per_cell = function(plot_pdf = T) {
  
  fig_dir = "figs/paper_figs/fig_s1/fig_s1a"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  n_umi = colSums(mat@mat)
  
  if(plot_pdf) {
    pdf("figs/paper_figs/fig_s1/fig_s1a/n_umi_per_cell.pdf",useDingbats = F)
  } else {
    png("figs/paper_figs/fig_s1/fig_s1/fig_s1a/n_umi_per_cell.png")
  }
  hist(n_umi,40,main = "UMI per cell",col = "grey80",xlab = "#UMI")
  dev.off()
  
}


fig_s1_a_number_of_cells_per_embryo = function() {
  
  fig_dir = "figs/paper_figs/fig_s1/fig_s1a"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  
  n_cls_per_emb = table(mat@cell_metadata[names(mc@mc),"transcriptional_rank"])
  
  pdf(sprintf("%s/number_of_cells_per_emb_gray.pdf",fig_dir),useDingbats = F)
  barplot(n_cls_per_emb,log = "y",col = "grey80")
  dev.off()
}

fig_s1_a_number_of_cls_per_age_group = function() {
  
  fig_dir = "figs/paper_figs/fig_s1/fig_s1a"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  
  emb_age = unique(mat@cell_metadata[names(mc@mc),c("transcriptional_rank","age_group","developmental_time")])
  
  time_age_group = tapply(emb_age$developmental_time,emb_age$age_group,median)
  
  n_cls_per_age_group = table(mat@cell_metadata[names(mc@mc),"age_group"])
  
  
  # estimated number of cells per age group
  emb_counts = read.table(file = "data/external_data/counted_embryos/wt10.cell_counts.txt",sep = "\t",stringsAsFactors = F)
  emb_time = unique(mat@cell_metadata[names(mc@mc),c("embryo","developmental_time")])
  
  emb_to_time = emb_time$developmental_time
  names(emb_to_time) = emb_time$embryo
  fit_count = lm(log2(emb_counts$cell_count) ~ emb_to_time[emb_counts$embryo])
  
  a = fit_count$coefficients[1]
  b = fit_count$coefficients[2]
  
  
  log_n_est_per_emb = a + b*time_age_group
  
  n_est_per_age_group = 2^log_n_est_per_emb
  
  mat_n_cls = cbind(n_cls_per_age_group,n_est_per_age_group)
  
  pdf(sprintf("%s/number_of_cells_per_age_group.pdf",fig_dir),useDingbats = F)
  barplot(t(mat_n_cls),beside = T,col = c("grey30","grey80"))
  legend(x = "topleft",legend = c("estimated in embryo from \n this age group","sampled in age group"),pch = 15, col = c("grey80","grey30"))
  dev.off()
}


fig_s1c = function() {
  
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  emb_ranks = unique(mat@cell_metadata[names(mc@mc),c("embryo","transcriptional_rank","morphology_rank","mouse_type")])
  emb_ranks$color = as.character(ifelse(emb_ranks$mouse_type == "ICR","cornflowerblue","black"))
  
  
  dir_name = "figs/paper_figs/fig_s1"
  
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  cex.lab = 2
  cex.axis = 2
  cex.main = 2
  margins = c(5,5,5,5)
  cex = 3
  
  
  xlabel = "Transcriptional rank"
  ylabel = "Morphological rank"
  main_label = ""
  #xlabel = "intrinsic rank"
  #ylabel = "atlas rank"
  #main_label = "Atlas vs Intrinsic Rank"
  
  pdf(sprintf("%s/fig_s1c.pdf",dir_name),useDingbats = F)
  #png(sprintf("%s/atlas_rank_vs_intrinsic_rank.png",dir_name))
  par(mar = margins)
  plot(x = emb_ranks$transcriptional_rank, y = emb_ranks$morphology_rank,pch = 19,cex.lab = cex.lab,cex.axis = cex.axis,
       xlab = xlabel,ylab = ylabel,col = emb_ranks$color,
       main = main_label,cex.main = cex.main,cex = cex)
  legend(x = "topleft",legend = c("ICR","C57BL/6"),pch = 19,col = c("cornflowerblue","black"),cex = 2)
  dev.off()
  
}

fig_s1e = function(plot_pdf = T) {
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  emb_ranks = unique(mat@cell_metadata[names(mc@mc),c("embryo","transcriptional_rank","developmental_time","age_group")])
  
  age_group_cols = c(RColorBrewer::brewer.pal(n = 12,name = "Paired"),"cornflowerblue")
  
  emb_ranks$color = age_group_cols[emb_ranks$age_group]
  
  dir_name = "figs/paper_figs/fig_s1"
  
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  cex.lab = 2
  cex.axis = 2
  cex.main = 2
  margins = c(5,6,5,4)
  cex = 2
  
  
  xlabel = "Transcriptional rank"
  ylabel = "Developmental time"
  main_label = ""
  #xlabel = "intrinsic rank"
  #ylabel = "atlas rank"
  #main_label = "Atlas vs Intrinsic Rank"
  if (plot_pdf) {
    pdf(sprintf("%s/fig_s1e.pdf",dir_name),useDingbats = F)
    #png(sprintf("%s/atlas_rank_vs_intrinsic_rank.png",dir_name))
    par(mar = margins)
    plot(x = emb_ranks$transcriptional_rank, y = emb_ranks$developmental_time,pch = 19,cex.lab = cex.lab,cex.axis = cex.axis,
         xlab = xlabel,ylab = ylabel,col = emb_ranks$color,
         main = main_label,cex.main = cex.main,cex = cex)
    legend(x = "topleft",legend = c(1:13),title = "Age group",pch = 19,col = age_group_cols, cex = 1)
    dev.off()
    
  } else {
    png(sprintf("%s/fig_s1e.png",dir_name))
    #png(sprintf("%s/atlas_rank_vs_intrinsic_rank.png",dir_name))
    par(mar = margins)
    plot(x = emb_ranks$transcriptional_rank, y = emb_ranks$developmental_time,pch = 19,cex.lab = cex.lab,cex.axis = cex.axis,
         xlab = xlabel,ylab = ylabel,col = emb_ranks$color,
         main = main_label,cex.main = cex.main,cex = cex)
    legend(x = "topleft",legend = c(1:13),title = "Age group",pch = 19,col = age_group_cols, cex = 1)
    dev.off()
    
  }

  
}



fig_s1f = function(plot_pdf = T) {
  mat_id = "sing_emb_wt10"
  mc_id = "sing_emb_wt10_recolored"
  dir_name = "figs/paper_figs/fig_s1"
  m_0 = 0.01
  s_0 = 0.005
  
  m = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  mc2d = scdb_mc2d("sing_emb_wt10_recolored")
  
  col_to_rank = c(1:nrow(mc@color_key))
  names(col_to_rank) = mc@color_key$color
  #mc_exe = c(which(mc@mc_fp["Ttr",] > 4),379)
  mc_exe = which(mc@colors %in% c("#7F6874","#F6BFCB","#EF5A9D","#F397C0","#0F4A9C"))
  
  
  cls_exe_all = names(mc@mc)[mc@mc %in% mc_exe]
  
  cls_exe = sample(x = cls_exe_all,500)
  
  
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
  
  max_s_score = quantile(s_score,0.9995)
  max_m_score = quantile(m_score,0.9995)
  
  f = s_score < max_s_score & m_score < max_m_score
  
  p_coldens =densCols(x = s_score,y = m_score,colramp = colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"),bias = 1))
  
  if(plot_pdf) {
    pdf(sprintf("%s/fig_s1f.pdf", dir_name), w=6, h=6,useDingbats = F)
    plot(s_score[f], m_score[f], pch=19, main = "S phase vs M phase UMIs",cex=1,
         xlab = "S phase score",ylab = "M phase score",
         col = p_coldens[f])
    points(s_score[cls_exe], m_score[cls_exe], pch=19, cex=0.8, col="black")
    abline(a = m_0,b = -m_0/s_0)
    dev.off()
  } else {
    png(sprintf("%s/fig_s1f.png", dir_name), w=2000, h=2000)
    plot(s_score[f], m_score[f], pch=19, main = "",cex=3,
         xlab = "",ylab = "",cex.axis = 2.5,
         col = p_coldens[f])
    points(s_score[cls_exe], m_score[cls_exe], pch=19, cex=2.4, col="black")
    abline(a = m_0,b = -m_0/s_0)
    dev.off()
  }
  
}



fig_s1g = function(plot_pdf = T) {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  legc = log2(mc@e_gc + 1e-5)
  
  
  # first separation embryonic endoderm (including node/notochord) from meso/-ectoderm 
  x1 = -16
  y1 = -12
  x2 = -12
  y2 = -16
  
  b_emb_endo = (y2 - y1)/(x2 - x1)
  a_emb_endo = (y1*x2 - y2*x1)/(x2 - x1)
  
  if(plot_pdf) {
    pdf(file = "figs/paper_figs/fig_s1/fig_s1g_Foxa2_Foxa1.pdf",useDingbats = F)
  } else {
    png(filename = "figs/paper_figs/fig_s1/fig_s1g_Foxa2_Foxa1.png")
  }
  par(mar = c(4,5,2,2))
  plot(x = legc["Foxa2",], y= legc["Foxa1",],pch = 19,col = mc@colors,cex = 2,xlab = "Foxa2",
       ylab = "Foxa1",cex.lab = 1)
  abline(a = a_emb_endo, b = b_emb_endo,lty = "dashed")
  dev.off()
  
  # second separation extraembryonic from embryonic endoderm
  x1 = -8.4
  y1 = -14
  x2 = -11
  y2 = -8.4
  
  b_exe_endo = (y2 - y1)/(x2 - x1)
  a_exe_endo = (y1*x2 - y2*x1)/(x2 - x1)
  
  if(plot_pdf) {
    pdf(file = "figs/paper_figs/fig_s1/fig_s1g_Apoe_Ttr.pdf",useDingbats = F)
  } else {
    png(filename = "figs/paper_figs/fig_s1/fig_s1g_Apoe_Ttr.png")
  }
  par(mar = c(4,5,2,1))
  plot(x = legc["Apoe",], y= legc["Ttr",],pch = 19,col = mc@colors,cex = 2,xlab = "Apoe",
       ylab = "Ttr",cex.lab = 1)
  abline(a = a_exe_endo, b = b_exe_endo,lty = "dashed")
  dev.off()
  
}



fig_s1h = function() {
  mat_id = "sing_emb_wt10"
  mc_id = "sing_emb_wt10_recolored"
  dir_name = "figs/paper_figs/fig_s1"
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
  
  s_genes = intersect(rownames(mc@mc_fp), s_genes)
  m_genes = intersect(rownames(mc@mc_fp), m_genes)
  
  tot  = colSums(m@mat)
  s_tot = colSums(m@mat[s_genes,])
  m_tot = colSums(m@mat[m_genes,])
  
  s_score = s_tot/tot
  m_score = m_tot/tot
  
  # boxplot per cell type of the cell cycle score
  
  cc_score_all = m_score + s_score
  
  all_cls = intersect(colnames(m@mat),names(mc@mc))
  
  cc_score = split(cc_score_all[all_cls],f = mc@colors[mc@mc[all_cls]])
  cc_score = cc_score[order(col_to_rank[names(cc_score)])]
  
  pdf(sprintf("%s/fig_s1h.pdf",dir_name),useDingbats = F)
  boxplot(cc_score,pch = 19,cex = 0.5,col = names(cc_score),xaxt = 'n')
  dev.off()
  
}



fig_s1i = function(plot_pdf = T) {
  
  dir_name = "figs/paper_figs/fig_s1"
  fn = sprintf("%s/fig_s1i_fraction_of_exe_endo_cells.png",dir_name)
  fn2 = sprintf("%s/fig_s1i_fraction_of_endoderm_cells.png",dir_name)
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  
  mc_exe = which(mc@colors %in% c("#7F6874","#F6BFCB"))
  
  mc_ag = table(mc@mc,mat@cell_metadata[names(mc@mc),"age_group"])
  mc_ag_n = t(t(mc_ag)/colSums(mc_ag))
  
  fr_exe_endo = colSums(mc_ag_n[mc_exe,])
  
  fit_y_exe_endo = 0.2*(0.83^c(0:12))
  
  if(plot_pdf) {
    pdf(file= gsub(".png",".pdf",fn),w = 8, h = 4,useDingbats = F)
  } else {
    png(filename = fn,w = 800, h = 400)
  }
  plot(x = c(1:13),y = fr_exe_endo,pch = 19,log = "y",xlab = "",ylab = "",cex = 4,cex.axis = 2)
  lines(x = c(1:13),fit_y_exe_endo)
  dev.off()
  
  
  mc_endo = which(mc@colors %in% c("#0F4A9C","#F397C0","#EF5A9D"))
  
  fr_endo = pmax(colSums(mc_ag_n[mc_endo,]),0.02)
  fit_y_endo = 0.10*(0.88^c(0:6))
  
  if(plot_pdf) {
    pdf(file = gsub(".png",".pdf",fn2),w = 8, h = 4,useDingbats = F)
  } else {
    png(filename = fn2,w = 800, h = 400)
  }
  plot(x = c(1:13),y = fr_endo,pch = 19,log = "y",xlab = "",ylab = "",cex = 4,cex.axis = 2)
  lines(x = c(7:13),fit_y_endo)
  dev.off()
  
  
  
}





