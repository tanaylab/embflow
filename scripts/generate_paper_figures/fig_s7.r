library("Matrix")
library("gridExtra")
source("scripts/foxc12/differential_expression_analysis.r")
source("scripts/foxc12/foxc_chim_timing.r")
source("scripts/foxc12/control_chim_timing.r")
source("scripts/foxc12/generate_chimera_tetraploid_data_analysis.r")
source("scripts/foxc12/definition_cell_types_endoderm_ectoderm_embryonic_mesoderm.r")

gen_fig_s7_plots = function() {
  
  # all the necessary data for this figure is saved in data/chimera_tetraploid_analysis/
  # if not available, please run the two functions below
  # (defined in scripts/foxc12/generate_chimera_tetraploid_data_analysis.r)
  # that creates all the necessary data for the Foxc DKO and Control tetraploid embryos
  #
  # foxc_tetraploid_generate_time_and_cell_type_annotation()
  # control_tetraploid_generate_time_and_cell_type_annotation()
  
  
  if(!dir.exists("figs/paper_figs/fig_s7")) {
    dir.create("figs/paper_figs/fig_s7")
  }
  
  # Figure S7C, like Figure 6E, relies additionally on the list of differentially expressed genes
  # saved in data/chimera_tetraploid_analysis/fig_6e_differentially_expressed_genes.txt
  # if not available run select_differentially_expressed_genes() to generate this list of genes
  # save it to data/chimera_tetraploid_analysis/fig_6e_differentially_expressed_genes.txt as tab-separated table
  # function select_differentially_expressed_genes() is defined in scripts/foxc12/differential_expression_analysis.r
  #
  # genes_f = select_differentially_expressed_genes()
  
  # Figure S7A: see Figure 6C
  
  fig_s7b()
  fig_s7c()
  fig_s7e()
  fig_s7f()
  
}

fig_s7b = function() {
  
  # Endo-/Ectoderm
  endo_colors = endo_ct_colors()
  ecto_colors = ectoderm_ct_colors()
  included_colors = c(endo_colors,ecto_colors)
  ks_ko_endo_ecto = foxc_chimera_plot_cumulative_distribution_ko_host(mat_nm = "foxc_chim_wt10",tag = "endo_ecto")
  ks_control_endo_ecto = control_chimera_plot_cumulative_distribution_control_host(mat_nm = "control_chim_wt10",tag = "endo_ecto")
  
  
  # Embryonic mesoderm
  included_colors = embryonic_meso_ct_colors()
  #foxc_chim_timing(tag = "emb_meso",included_colors = included_colors)
  ks_ko_emb_meso = foxc_chimera_plot_cumulative_distribution_ko_host(mat_nm = "foxc_chim_wt10",tag = "emb_meso")
  ks_control_emb_meso = control_chimera_plot_cumulative_distribution_control_host(mat_nm = "control_chim_wt10",tag = "emb_meso")
  
  
  y_endo = c(ks_ko_endo_ecto$d_ks,ks_control_endo_ecto$d_ks)
  names(y_endo) = c(ks_ko_endo_ecto$chimera,ks_control_endo_ecto$chimera)
  
  y_emb_meso = c(ks_ko_emb_meso$d_ks,ks_control_emb_meso$d_ks)
  names(y_emb_meso) = c(ks_ko_emb_meso$chimera,ks_control_emb_meso$chimera)
  
  
  chim_time_ko = read.table("data/chimera_tetraploid_analysis/foxc_chim_wt10/time_match/time_match_summary.txt",sep = "\t",stringsAsFactors = F,h = T)
  chim_time_control = read.table("data/chimera_tetraploid_analysis/control_chim_wt10/time_match/time_match_summary.txt",sep = "\t",stringsAsFactors = F,h = T)
  chim_time_ko = chim_time_ko[order(chim_time_ko$best_rank_host),]
  chim_time_control = chim_time_control[order(chim_time_control$best_rank_host),]
  
  rownames(chim_time_ko) = chim_time_ko$embryo
  rownames(chim_time_control) = chim_time_control$embryo
  
  emb_col = c(rep("#83c26d",length(ks_ko_emb_meso$delta_median)),rep("cornflowerblue",length(ks_control_emb_meso$delta_median)))
  
  
  pdf("figs/paper_figs/fig_s7/fig_s7b.pdf",useDingbats = F)
 
  plot(x =  y_endo[c(chim_time_ko$embryo,chim_time_control$embryo)],
       y =  y_emb_meso[c(chim_time_ko$embryo,chim_time_control$embryo)],
       pch = 19,cex = 2,
       col = emb_col,
       xlab = "D(KS) Enoderm/Ectoderm",
       ylab = "D(KS) Embryonic mesoderm",     
       xlim = c(0,0.35),
       ylim = c(0,0.35))
  abline(a = 0,b = 1,lty = "dashed")
  points(x =  y_endo[c(chim_time_ko$embryo,chim_time_control$embryo)],
         y =  y_emb_meso[c(chim_time_ko$embryo,chim_time_control$embryo)],
         pch = 19,cex = 3,
         col = emb_col)
  dev.off()
  
  
  
  
}

fig_s7c = function() {
  
  reg = 5e-5
  
  included_colors = embryonic_meso_ct_colors()
  
  genes_f = read.table("data/chimera_tetraploid_analysis/fig_6e_differentially_expressed_genes.txt",sep = "\t",stringsAsFactors = F)$x
  fig_dir = "figs/paper_figs/fig_s7/fig_s7c"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mat_nm = "control_chim_wt10"
  mat = scdb_mat(mat_nm)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  
  load(file = sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
  query_cls_col = cmp_annot$query_cls_col
  
  df_chim = read.table(sprintf("data/chimera_tetraploid_analysis/%s/time_match/time_match_summary.txt",mat_nm),sep = "\t",h = T,stringsAsFactors = F)
  best_wt_rank = df_chim$best_rank_host
  best_wt_rank = round(best_wt_rank)
  names(best_wt_rank) = df_chim$embryo
  
  control_cls = names(query_cls_col)[( mat@cell_metadata[names(query_cls_col),"cell_type"] == "control" ) & ( query_cls_col %in% included_colors )]
  
  cmp_emb = cmp_differential_expression_to_wt_per_emb(mat_query = mat,
                                                      query_cls = control_cls,
                                                      chim_best_wt_rank = best_wt_rank,
                                                      included_colors = included_colors)
  
  lfp_emb = log2(cmp_emb$egc_query + reg) - log2(cmp_emb$egc_wt + reg)
  
  ct_to_col = mc_wt@color_key$color
  names(ct_to_col) = mc_wt@color_key$group
  
  ct_order = c("Caudal mesoderm","Paraxial mesoderm",
               "Rostral mesoderm","Cardiac mesoderm","Lateral & intermediate mesoderm",
               "Amnion/Chorion","ExE mesoderm","Allantois")
  
  included_colors_ct = ct_to_col[ct_order]
  
  control_cls_ct = names(query_cls_col)[( mat@cell_metadata[names(query_cls_col),"cell_type"] == "control" ) & ( query_cls_col %in% included_colors_ct )]
  
  cmp_ct = cmp_differential_expression_to_wt_per_ct(mat_query = mat,
                                                    query_cls_col = query_cls_col[control_cls_ct],
                                                    chim_best_wt_rank = best_wt_rank,
                                                    included_colors = included_colors_ct)
  
  lfp_ct = log2(cmp_ct$egc_query + reg) - log2(cmp_ct$egc_wt + reg)
  
  lfp_ct = lfp_ct[genes_f,ct_to_col[ct_order]]
  lfp_emb = lfp_emb[genes_f,]
  
  plot_heatmap_diff_expression_per_emb_and_ct(lfp_emb = lfp_emb,lfp_ct = lfp_ct,fig_dir = fig_dir,plot_pdf = T,tag = "chim_control")
  
  
}

fig_s7d = function() {
  
  mat_nm1 = "foxc_tetra_wt10"
  mat_nm2 = "control_tetra_wt10"
  
  foxc_col = "#83c26d"
  control_col = "cornflowerblue"
  
  fig_scale = 1.4
  lwd = 8
  xlim_min = 1
  xlim_max = 153
  
  load("data/chimera_tetraploid_analysis/foxc_tetra_wt10/time_match/time_dist_ko.Rda")
  load("data/chimera_tetraploid_analysis/control_tetra_wt10/time_match/time_dist_control.Rda")
  
  ko_dist_all = time_dist_ko$ko
  control_dist_all = time_dist_control$control
  
  pdf("figs/paper_figs/fig_s7/fig_s7d.pdf",w = 5*fig_scale,h = 4*fig_scale,useDingbats = F)
  ko_dens = cumsum(ko_dist_all[1,])/sum(ko_dist_all[1,])
  plot(c(xlim_min:xlim_max),ko_dens[xlim_min:xlim_max],type = "l",
       lwd = lwd,col = foxc_col,ylim = c(0,1),xlab = "",ylab = "",
       main = "Time Distribution 4N Embryos")
  
  for (i in 2:nrow(ko_dist_all)) {
    ko_dens = cumsum(ko_dist_all[i,])/sum(ko_dist_all[i,])
    lines(c(xlim_min:xlim_max),ko_dens[xlim_min:xlim_max],type = "l",lwd = lwd,col = foxc_col)
  }
  for (i in 1:nrow(control_dist_all)) {
    control_dens = cumsum(control_dist_all[i,])/sum(control_dist_all[i,])
    lines(c(xlim_min:xlim_max),control_dens[xlim_min:xlim_max],type = "l",lwd = lwd,col = control_col)
  }
  dev.off()
  
  
}



fig_s7e = function() {
  
  fig_dir= "figs/paper_figs/fig_s7/fig_s7e"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  barplot_ct_frequency(mat_nm = "foxc_tetra_wt10",fig_dir = fig_dir,plot_pdf = T)
  barplot_ct_frequency(mat_nm = "control_tetra_wt10",fig_dir = fig_dir,plot_pdf = T)

}


barplot_ct_frequency = function(mat_nm,fig_dir,plot_pdf = T) {
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  df_chim = read.table(sprintf("data/chimera_tetraploid_analysis/%s/time_match/time_match_summary.txt",mat_nm),sep = "\t",stringsAsFactors = F,h= T)
  rownames(df_chim) = df_chim$embryo
  mat = scdb_mat(mat_nm)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  col_to_rank = c(1:nrow(mc_wt@color_key))
  names(col_to_rank) = mc_wt@color_key$color
  
  excluded_colors = c("#F6BFCB","#7F6874") 
  included_colors = setdiff(unique(mc_wt@color_key$color),excluded_colors)
  
  
  query_embryos = df_chim$embryo
  
  age_field = c("best_rank_host","best_rank_host","best_rank_ko","best_rank_control")
  names(age_field) = c("foxc_chim_wt10","control_chim_wt10","foxc_tetra_wt10","control_tetra_wt10")
  
  query_embryos = query_embryos[order(df_chim[query_embryos,age_field[mat_nm]])]
  
  load(file = sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
  
  query_cls_col = cmp_annot$query_cls_col
  query_cls = names(query_cls_col)[!(query_cls_col %in% excluded_colors)]
  query_cls = query_cls[mat@cell_metadata[query_cls,"embryo"] %in% query_embryos]
  
  tmp = matrix(0,nrow = length(query_embryos),ncol = length(included_colors))
  rownames(tmp) = query_embryos
  colnames(tmp) = included_colors
  
  query_cls_f = query_cls[( mat@cell_metadata[query_cls,"cell_type"] %in% c("KO","control") ) & ( mat@cell_metadata[query_cls,"embryo"] %in% query_embryos)]
  
  query_vs_ct_tmp = table(mat@cell_metadata[query_cls_f,"embryo"],query_cls_col[query_cls_f])
  query_vs_ct = tmp
  query_vs_ct[rownames(query_vs_ct_tmp),colnames(query_vs_ct_tmp)] = query_vs_ct_tmp
  query_vs_ct = query_vs_ct[query_embryos,order(col_to_rank[colnames(query_vs_ct)])] 
  
  n_cls_ko = rowSums(query_vs_ct)
  query_vs_ct_n = query_vs_ct/rowSums(query_vs_ct)
  
  if(plot_pdf) {
    pdf(sprintf("%s/%s_ct_freq.pdf",fig_dir,mat_nm))
    barplot(t(query_vs_ct_n),col = colnames(query_vs_ct_n),las = 2,axisnames = F)
    dev.off()
  } else {
    png(sprintf("%s/%s_ct_freq.png",fig_dir,mat_nm),w = 550,h= 550)
    barplot(t(query_vs_ct_n),col = colnames(query_vs_ct_n),las = 2,axisnames = F)
    dev.off()
  }
  
  
}


fig_s7f = function() {
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  
  p1 = chimera_ko_vs_host_expr(ct = mc_wt@color_key$color[1])
  p2 = chimera_ko_vs_host_expr(ct = mc_wt@color_key$color[10])
  p3 = tetraploid_ko_vs_wt_expr(ct = mc_wt@color_key$color[1])
  p4 = tetraploid_ko_vs_wt_expr(ct = mc_wt@color_key$color[10])
  
  pdf("figs/paper_figs/fig_s7/fig_s7f.pdf",h = 10,w = 10,useDingbats = F)
  grid.arrange(grobs = list(p1,p2,p3,p4),nrow = 2,ncol = 2,heights = c(1,1),widths = c(1,1))
  dev.off()
  
}


chimera_ko_vs_host_expr = function(ct,lfp_thr = 1.5) {
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  mat_wt = scdb_mat("sing_emb_wt10")
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  
  reg = 1e-4
  max_non_meso_lfp = 1
  
  foxc_chim_highlighted_genes = read.table("data/chimera_tetraploid_analysis/fig_6e_differentially_expressed_genes.txt",sep = "\t",stringsAsFactors = F)$x
  
  bad_genes = read.table(file = "data/external_data/sing_emb_wt10.bad_genes.txt",sep = "\t",stringsAsFactors = F)$x
  igf2_correlated_genes = c("Igf2","AK145379;H19","Polg","Slc25a4","Peg10","Igf2as","AK086477;Sh3glb1")
  
  feat_genes = setdiff(rownames(mc_wt@e_gc),bad_genes)
  
  # next FoxC chimera
  
  mat_nm_chim = "foxc_chim_wt10"
  mat_chim= scdb_mat(mat_nm_chim)
  
  
  load(file = sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm_chim))
  ko_cls = colnames(mat_chim@mat)[(mat_chim@cell_metadata[colnames(mat_chim@mat),"cell_type"] == "KO") & (cmp_annot$query_cls_col[colnames(mat_chim@mat)] %in% ct)]
  host_cls = colnames(mat_chim@mat)[(mat_chim@cell_metadata[colnames(mat_chim@mat),"cell_type"] == "host") & (cmp_annot$query_cls_col[colnames(mat_chim@mat)] %in% ct)]
  
  egc_ko = rowSums(mat_chim@mat[feat_genes,ko_cls])
  egc_host = rowSums(mat_chim@mat[feat_genes,host_cls])
  egc_ko = egc_ko/sum(egc_ko)
  egc_host = egc_host/sum(egc_host)
  
  lfp_chim = log2(egc_ko + reg) - log2(egc_host + reg)
  
  egc_ko = log2(egc_ko + reg)
  egc_host =  log2(egc_host + reg)
  
  gene_highlight = names(lfp_chim) %in% foxc_chim_highlighted_genes
  
  df_col = ifelse(gene_highlight,"purple","gray")
  
  gene_diff = abs(lfp_chim) > lfp_thr
  
  gene_highlight = gene_highlight | gene_diff
  
  df_col[gene_diff] = "red"
  df_col[names(lfp_chim) %in% foxc_chim_highlighted_genes] = "purple"

  df = data.frame(ko =egc_ko, host=egc_host, gene = names(lfp_chim), col=df_col,highlight = gene_highlight)
  
  min_egc = min(egc_host,egc_ko)
  max_egc = max(egc_host,egc_ko)
  
  df$gene = gsub("AK145379;H19","H19",x = df$gene)
  
  p = ggplot(df, aes(x = ko,y = host)) +
    geom_point(col=df$col) +
    ggtitle(paste0(paste(col_to_ct[ct],collapse = " ")," 2N vs host")) +
    geom_text_repel(data = filter(df, highlight),
                    aes(label=gene)) + 
    geom_point(data = filter(df,col == "purple"),col = "purple") +
    xlim(c(min_egc,max_egc)) +
    ylim(c(min_egc,max_egc)) +
    xlab("Foxc DKO (2N)") + 
    ylab("Host")
  
  ct_name = ifelse(col_to_ct[ct] == "Amnion/Chorion","Amnion_Chorion",col_to_ct[ct])
  
  if(length(ct_name) > 1) {
    ct_name = paste(ct_name,collapse = "_")
  }
  
  
  return(p)
}


tetraploid_ko_vs_wt_expr = function(ct,lfp_thr = 1.5) {
  
  
  mat_nm = "foxc_tetra_wt10"
  mat = scdb_mat(mat_nm)
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  mat_wt = scdb_mat("sing_emb_wt10")
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  
  load(file = sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
  cls_a = colnames(mat@mat)[colnames(mat@mat) %in% names(cmp_annot$query_cls_col)]
  ko_cls = cls_a[(mat@cell_metadata[cls_a,"cell_type"] == "KO") & (cmp_annot$query_cls_col[cls_a] %in% ct)]
  
  wt_cls = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% ct]
  
  reg = 1e-4
  max_non_meso_lfp = 1
  
  foxc_chim_highlighted_genes = read.table("data/chimera_tetraploid_analysis/fig_6e_differentially_expressed_genes.txt",sep = "\t",stringsAsFactors = F)$x
  
  bad_genes = read.table(file = "data/external_data/sing_emb_wt10.bad_genes.txt",sep = "\t",stringsAsFactors = F)$x
  igf2_correlated_genes = c("Igf2","AK145379;H19","Polg","Slc25a4","Peg10","Igf2as","AK086477;Sh3glb1")
  
  feat_genes = setdiff(rownames(mc_wt@e_gc),bad_genes)
  
  egc_query_tetra = rowSums(mat@mat[feat_genes,ko_cls])/sum(mat@mat[feat_genes,ko_cls])
  egc_ref_tetra = rowSums(mat_wt@mat[feat_genes,wt_cls])/sum(mat_wt@mat[feat_genes,wt_cls])
  
  lfp_tetra = log2(egc_query_tetra + reg) - log2(egc_ref_tetra + reg)
  
  egc_query_tetra = log2(egc_query_tetra + reg)
  egc_ref_tetra = log2(egc_ref_tetra + reg)
  
  gene_highlight = names(lfp_tetra) %in% foxc_chim_highlighted_genes
  
  df_col = ifelse(gene_highlight,"purple","gray")
  
  gene_diff = abs(lfp_tetra) > lfp_thr
  
  gene_highlight = gene_highlight | gene_diff
  
  df_col[gene_diff] = "red"
  df_col[names(lfp_tetra) %in% foxc_chim_highlighted_genes] = "purple"
  
  df = data.frame(tetraploid=egc_query_tetra, wt=egc_ref_tetra, gene = names(lfp_tetra), col=df_col,highlight = gene_highlight)
  
  df$gene = gsub("AK145379;H19","H19",x = df$gene)
  
  min_egc = min(egc_query_tetra,egc_ref_tetra)
  max_egc = max(egc_query_tetra,egc_ref_tetra)
  
  
  p = ggplot(df, aes(x = tetraploid, y = wt)) +
    geom_point(col=df$col) +
    ggtitle(paste(paste(col_to_ct[ct],collapse = " "),"4N vs WT"," ")) +
    geom_text_repel(data = filter(df, highlight),
                    aes(label=gene)) + 
    geom_point(data = filter(df,col == "purple"),col = "purple") +
    xlim(c(min_egc,max_egc)) +
    ylim(c(min_egc,max_egc)) +
    xlab("Foxc DKO (4N)") + 
    ylab("WT Expression")
  
  ct_name = ifelse(col_to_ct[ct] == "Amnion/Chorion","Amnion_Chorion",col_to_ct[ct])
  
  if(length(ct_name) > 1) {
    ct_name = paste(ct_name,collapse = "_")
  }
  
  return(p) 
}


