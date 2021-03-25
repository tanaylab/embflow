library("Matrix")
source("scripts/generate_paper_figures/plot_network.r")
source("scripts/foxc12/preprocessing/foxc_control_gating.r")

gen_fig_s6_plots = function() {
  
  if(!dir.exists("figs/paper_figs/fig_s6")) {
    dir.create("figs/paper_figs/fig_s6")
  }
  
  fig_s6a()
  fig_s6b()
  fig_s6c()
}

fig_s6b = function() {
  
  fig_dir = "figs/paper_figs/fig_s6/fig_s6b"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  # from scripts/foxc12/foxc_control_gating
  # this script was used to gate the single cells and assign each cell to one of
  # the groups KO,control,host,wt or unclear
  # this information was added to the metadata entry of the metacell mat object.
  # See end of the function gating_of_foxc_control()
  df_gating = gating_of_foxc_control(fig_dir)
  
}

fig_s6a = function() {
  
  mct_id = "sing_emb_wt10"
  mat_id = "sing_emb_wt10"
  genes = c("Foxc1","Foxc2")
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  mat_n = t(t(mat@mat[genes,names(mc@mc)])/colSums(mat@mat[,names(mc@mc)]))
  colors_ordered = mc@color_key$color
  
  cls_spl = split(names(mc@mc),mat@cell_metadata[names(mc@mc),"age_group"])
  
  reg = 3e-5
  
  legc_ls = lapply(genes,function(gene) {
    e_gc_mat = sapply(cls_spl,function(cls) {
      e_gc_tmp = tapply(mat_n[gene,cls],mc@mc[cls],mean)
      e_gc_t = rep(0,ncol(mc@e_gc))
      e_gc_t[as.numeric(names(e_gc_tmp))] = e_gc_tmp
      return(e_gc_t)
    })
    legc = log2(reg + e_gc_mat)
    return(legc)
  })
  
  
  mm_mctnetwork_plot_net(mct_id = mct_id,fn = "figs/paper_figs/fig_s6/fig_s6a_foxc1_over_flows.png",colors_ordered = colors_ordered,mc_t_score = legc_ls[[1]],dy_ext = 0,dx_back = 0,w = 2200,h = 1200)
  mm_mctnetwork_plot_net(mct_id = mct_id,fn = "figs/paper_figs/fig_s6/fig_s6a_foxc2_over_flows.png",colors_ordered = colors_ordered,mc_t_score = legc_ls[[2]],dy_ext = 0,dx_back = 0,w = 2200,h = 1200)
}

fig_s6a_color_scale = function() {
  
  cols = colorRampPalette(c("lightgray", "gray", "darkgray", "lightpink", "pink", "red", "darkred"))(101)
  vals = seq(0,1,length.out = 101)
  show_vals_ind = c(1,51,101)
  
  pdf(file = "figs/paper_figs/fig_s6/fig_s6a_color_scale.pdf",useDingbats = F)
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')
  
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  dev.off()
  
}



fig_s6c = function(plot_pdf = T) {
  
  fig_dir = "figs/paper_figs/fig_s6/fig_s6c"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  genes_mm9 = read.table("data/external_data/gene_intervals_mm9.txt",sep = "\t",stringsAsFactors = F,header = T)
  mat=  scdb_mat("foxc_chim_wt10")
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  cgraph = scdb_cgraph("foxc_chim_wt10")
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  col_to_ct = c(col_to_ct,"Mixed")
  names(col_to_ct)[length(col_to_ct)] = "gray"
  included_colors = setdiff(unique(mc_wt@colors),c("#F6BFCB","#7F6874"))
  
  df_chim = read.table("data/chimera_tetraploid_analysis/foxc_chim_wt10/time_match/time_match_summary.txt",h= T,sep = "\t",stringsAsFactors = F)
  chim_embryos = df_chim$embryo
  rownames(df_chim) = df_chim$embryo
  chim_embryos = chim_embryos[order(df_chim[chim_embryos,"best_rank_host"])]
  ref_ranks = c()
  for (i in 1:length(chim_embryos)) {
    host_match = round(df_chim[chim_embryos[i],"best_rank_host"])
    ref_ranks = c(ref_ranks,c((host_match-2):(host_match+2)))
  }
  ref_ranks = unique(ref_ranks)
  
  load("data/chimera_tetraploid_analysis/foxc_chim_wt10/color_annotation/cmp_annot.Rda")
  query_cls_col = cmp_annot$query_cls_col
  
  all_cells = colnames(mat@mat)
  
  host_cls = all_cells[(mat@cell_metadata[all_cells,"cell_type"] == "host") & (mat@cell_metadata[all_cells,"embryo"] %in% chim_embryos)]
  host_cls = host_cls[query_cls_col[host_cls] %in% included_colors]
  
  fox_ko_cls =  all_cells[mat@cell_metadata[all_cells,"cell_type"] == "KO" & (mat@cell_metadata[all_cells,"embryo"] %in% chim_embryos)]
  fox_ko_cls = fox_ko_cls[query_cls_col[fox_ko_cls] %in% included_colors]                       
  
  wt10_cls = names(mc_wt@mc)[(mat@cell_metadata[names(mc_wt@mc),"transcriptional_rank"] %in% ref_ranks) & mc_wt@colors[mc_wt@mc] %in% included_colors]
  wt10_cls = intersect(wt10_cls,colnames(mat@mat))
  
  egc_ko = rowSums(mat@mat[,fox_ko_cls])/sum(mat@mat[,fox_ko_cls])
  egc_host = rowSums(mat@mat[,host_cls])/sum(mat@mat[,host_cls])
  egc_ref = rowSums(mat@mat[,wt10_cls])/sum(mat@mat[,wt10_cls])
  
  
  lfp_ko_host = (egc_ko + 1e-4)/(egc_host + 1e-4)
  lfp_ko_ref = (egc_ko + 1e-4)/(egc_ref + 1e-4)
  lfp_host_ref = (egc_host + 1e-4)/(egc_ref + 1e-4)
  lfp_ko_host = (egc_ko + 1e-5)/(egc_host + 1e-5)
  lfp_ko_ref = (egc_ko + 1e-5)/(egc_ref + 1e-5)
  lfp_host_ref = (egc_host + 1e-5)/(egc_ref + 1e-5)
  
  included_chromosomes = grep("chr",unique(genes_mm9$chrom),v = T)
  included_chromosomes = setdiff(included_chromosomes,"chrM")
  
  mm9_f = genes_mm9[genes_mm9$chrom %in% included_chromosomes,]  
  
  gene_to_chr =   unique(mm9_f[,c("chrom","gene_name")])
  
  gene_freq = table(gene_to_chr$gene_name)
  genes_f = names(gene_freq)[gene_freq == 1]
  
  mm9_f = mm9_f[mm9_f$gene_name %in% genes_f,]
  
  gene_pos = tapply(mm9_f$start,mm9_f$gene_name,min)
  
  gene_to_chr = gene_to_chr[gene_to_chr$gene_name %in% genes_f,]
  
  gene_ls = split(gene_to_chr$gene_name,gene_to_chr$chrom)
  chromosomes = paste0("chr",as.character(c(c(1:19),c("X","Y"))))
  gene_ls =gene_ls[chromosomes]
  
  ymin = 0
  ymax = 2.5
  abline_h = 1
  if(plot_pdf) {
    h_all = 10.5
    h_one = 0.45
    h_first = 0.8
    w_all = 15
  } else {
    
    h_all = 1050
    h_one = 45
    h_first = 80
    w_all = 1500
  }
  
  if(plot_pdf) {
    pdf(sprintf("%s/karyogramm_ko_vs_host.pdf",fig_dir),h = h_all, w = w_all,useDingbats = F)
  } else {
    png(sprintf("%s/karyogramm_ko_vs_host.png",fig_dir),h = h_all, w = w_all)
  }
  
  layout(mat = matrix(c(1:21),nrow = 21,ncol = 1),heights = c(h_first,rep(h_one,20)))
  par(mar= c(0.5,6,4,4))
  for(chrom in chromosomes) {
    genes = gene_ls[[chrom]]
    if(chrom == "chr1") {
      main_plot = "Chimera KO vs host"
    } else {
      main_plot = ""
    }
    plot(x = gene_pos[genes],y = pmax(pmin(lfp_ko_host[genes],ymax),ymin),pch =19,ylim = c(ymin,ymax),xlab = "",ylab = chrom,cex.lab = 2,main = main_plot,
         cex.main = 4)
    abline(h = abline_h)
    par(mar= c(0.5,6,0.5,4))
  }
  dev.off()
  
  for (chrom in c("chr1","chr7","chr8")) {
    genes = gene_ls[[chrom]]
    if(plot_pdf) {
      pdf(sprintf("%s/%s.pdf",fig_dir,chrom), w = 10,h = 3)
    } else {
      png(sprintf("%s/%s.png",fig_dir,chrom), w = 1000,h = 300)
    }
    
    main_plot = ""
    par(mar = c(4,6,2,2))
    plot(x = gene_pos[genes],y = pmax(pmin(lfp_ko_host[genes],ymax),ymin),pch =19,ylim = c(ymin,ymax),xlab = "",ylab = chrom,cex.lab = 2,main = main_plot,
         cex.main = 4)
    abline(h = abline_h)
    dev.off()
    
  }
  
}