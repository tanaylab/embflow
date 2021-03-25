library("tidyverse")
library("ggrepel")
library("Matrix")
library("zoo")
source("scripts/foxc12/definition_cell_types_endoderm_ectoderm_embryonic_mesoderm.r")



# the following function generates line graphs as shown in Figure 6F for a given group of genes gnms_f
plot_gene_chimera = function(df_plot,gnms_f,legc_wt,fig_dir,plot_pdf = F) {
  
  
  roll_width = 5
  
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  if(!dir.exists(sprintf("%s/genes",fig_dir))) {
    dir.create(sprintf("%s/genes",fig_dir))
  }
  
  xlim_min = 105
  
  for (gene in gnms_f) {
    print(gene)
    
    wt_ranks = intersect(c((xlim_min- roll_width):153),as.numeric(colnames(legc_wt)))
    
    legc_gene_wt = legc_wt[gene,as.character(wt_ranks)]
    legc_gene_wt = c(legc_gene_wt,legc_gene_wt[(length(legc_gene_wt) - roll_width + 1):length(legc_gene_wt)])
    
    mov_sd = rollapply(data = legc_gene_wt,width = 2*roll_width+1,sd)
    legc_gene_wt = rollmean(x = legc_gene_wt,k = 2*roll_width + 1)
    
    upper_freq = legc_gene_wt + mov_sd
    lower_freq = legc_gene_wt - mov_sd
    
    ylim_max = max(legc_gene_wt,df_plot[,gene])
    ylim_min = log2(1e-5)
    
    gene_nm = gsub(";","_",gene)
    if(plot_pdf) {
      
      pdf(sprintf("%s/genes/%s.pdf",fig_dir,gene_nm),w = 7,h = 5.6,useDingbats = F)
      par(mar = c(3,3,5,3))
      plot(x = df_plot[,"best_rank_host"],
           y = df_plot[,gene],
           col = df_plot$color,
           pch = df_plot$pch,
           main = gene,
           ylim = c(ylim_min,ylim_min + 1.1*(ylim_max - ylim_min)),
           cex = 4,cex.main = 4)
      polygon(x = c(wt_ranks[(roll_width+1):length(wt_ranks)],rev(wt_ranks[(roll_width+1):length(wt_ranks)])),y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
      lines(x = wt_ranks[(roll_width+1):length(wt_ranks)],y = legc_gene_wt,lwd = 6)
      points(x = df_plot[,"best_rank_host"],y = df_plot[,gene],col = df_plot$color,pch = df_plot$pch,cex = 4)
      
      dev.off()
      
    } else {
      
      png(sprintf("%s/genes/%s.png",fig_dir,gene_nm),w = 800,h = 600)
      par(mar = c(3,3,5,3))
      plot(x = as.numeric(colnames(legc_wt)),y = legc_wt[gene,],pch = 19,main = gene,ylim = c(ylim_min,ylim_min + 1.1*(ylim_max - ylim_min)),cex = 6,cex.main = 6)
      polygon(x = c(wt_ranks[(roll_width+1):length(wt_ranks)],rev(wt_ranks[(roll_width+1):length(wt_ranks)])),y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
    
      lines(x = wt_ranks[(roll_width+1):length(wt_ranks)],y = legc_gene_wt,lwd = 8)
      points(x = df_plot[,"best_rank_host"],y = df_plot[,gene],col = df_plot$color,pch = df_plot$pch,cex = 6)
      
      dev.off()
    }
    
    
    
  }
  
  
  
}


# computes differential expression between query and wt cells per embryo, 
# as shown in the left heatmap of Figure 6E and Figure S7C
cmp_differential_expression_to_wt_per_emb = function(mat_query,query_cls,chim_best_wt_rank,included_colors) {
  
  
  mat_wt = scdb_mat("sing_emb_wt10")
  mc_wt  = scdb_mc("sing_emb_wt10_recolored")
  
  wt10_cls = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% included_colors]
  
  query_cls_ls = split(query_cls,mat_query@cell_metadata[query_cls,"embryo"])
  
  egc_query = sapply(query_cls_ls,FUN = function(cls_emb) {
    return(rowSums(mat_query@mat[,cls_emb]))
  })
  tot_umi_query = colSums(egc_query)
  egc_query = t(t(egc_query)/colSums(egc_query))
  
  chim_best_wt_rank = chim_best_wt_rank[names(chim_best_wt_rank) %in% colnames(egc_query)]
  egc_query = egc_query[,names(chim_best_wt_rank)]
  
  egc_wt = sapply(1:length(chim_best_wt_rank),function(i) {
    
    wt_rank = chim_best_wt_rank[i]
    
    ref_ranks = c(max(wt_rank - 3,0):min(wt_rank+3,153))
    
    cls_f = wt10_cls[mat_wt@cell_metadata[wt10_cls,"transcriptional_rank"] %in% ref_ranks]
    
    egc = rowSums(mat_wt@mat[,cls_f])
    
    return(egc)
  })
  colnames(egc_wt) = names(chim_best_wt_rank)
  tot_umi_ref = colSums(egc_wt)
  egc_wt = t(t(egc_wt)/colSums(egc_wt))
  
  return(list(egc_query = egc_query,tot_umi_query = tot_umi_query,egc_wt = egc_wt))
}

# computes differential expression between query and wt cells per cell type, 
# as shown in the right heatmap of Figure 6E and Figure S7C
cmp_differential_expression_to_wt_per_ct = function(mat_query,query_cls_col,chim_best_wt_rank,included_colors) {
  
  query_cls = names(query_cls_col)
  
  mat_wt = scdb_mat("sing_emb_wt10")
  mc_wt  = scdb_mc("sing_emb_wt10_recolored")
  
  ref_ranks = c()
  for (wt_rank in chim_best_wt_rank) {
    ref_ranks = c(ref_ranks,c(max(wt_rank-2,0),min(wt_rank+2,153)))
  }
  ref_ranks = unique(ref_ranks)
  wt10_cls = names(mc_wt@mc)[( mc_wt@colors[mc_wt@mc] %in% included_colors ) & ( mat_wt@cell_metadata[names(mc_wt@mc),"transcriptional_rank"] %in% ref_ranks )]  
  
  
  query_cls_ls = split(query_cls,query_cls_col)
  
  
  egc_query = sapply(query_cls_ls,FUN = function(cls_ct) {
    return(rowSums(mat_query@mat[,cls_ct]))
  })
  tot_umi_query = colSums(egc_query)
  egc_query = t(t(egc_query)/colSums(egc_query))
  
  wt10_cls_ls = split(wt10_cls,mc_wt@colors[mc_wt@mc[wt10_cls]])
  egc_wt = sapply(wt10_cls_ls,FUN = function(cls_ct) {
    return(rowSums(mat_wt@mat[,cls_ct]))
  })
  egc_wt = t(t(egc_wt)/colSums(egc_wt))
  
  egc_wt = egc_wt[,colnames(egc_query)]
  
  
  return(list(egc_query = egc_query,tot_umi_query = tot_umi_query,egc_wt = egc_wt))
}


# plots heatmaps of Figure 6E and Figure S7C
plot_heatmap_diff_expression_per_emb_and_ct = function(lfp_emb,lfp_ct,fig_dir,plot_pdf = F,tag ="all") {

  # next cluster gnms_f
  lfp_clst = cbind(lfp_emb,lfp_ct)
  gg_hclst = hclust(as.dist(1 - tgs_cor(t(lfp_clst))))
  
  shades = colorRampPalette(rev(RColorBrewer::brewer.pal(9,"RdBu")))(100)
  breaks = c(seq(-2,0,length.out = 51)[-51],c(0),seq(0,2,length.out = 51)[-1])
  
  #lfp_emb = lfp_emb[gg_hclst$order,]
  #lfp_ct = lfp_ct[gg_hclst$order,]
  
  # plot the two heatmaps
  
  colnames(lfp_emb) = c(1:ncol(lfp_emb))
  
  if(plot_pdf) {
    
    pdf(sprintf("%s/diff_expr_per_embryo_and_cell_type_%s.pdf",fig_dir,tag),w = 10, h = 11,useDingbats = F)
    layout(matrix(c(1,0,2,3),nrow= 2,ncol = 2),widths = c(400,500),heights = c(900,80))
    
    
    par(mar= c(0,4,8,6))
    image(pmin(pmax(t(lfp_emb),-2),2),col = shades,breaks = breaks,axes = F)
    mtext(text = substr(rownames(lfp_emb),1,9),side = 4,at = seq(0,1,length.out = nrow(lfp_emb)),las = 2)
    
    mtext(text = colnames(lfp_emb),side = 1,at = seq(0,1,length.out = ncol(lfp_emb)),las = 2)
    
    
    par(mar= c(0,4,8,6))
    image(pmin(pmax(t(lfp_ct),-2),2),col = shades,breaks = breaks,axes = F)
    mtext(text = substr(rownames(lfp_emb),1,9),side = 4,at = seq(0,1,length.out = nrow(lfp_emb)),las = 2)
    
    par(mar = c(4,4,0,6))
    image(as.matrix(c(1:ncol(lfp_ct))),col = colnames(lfp_ct),axes = F)
    
    dev.off()
    
  } else {
    png(sprintf("%s/diff_expr_per_embryo_and_cell_type_%s.png",fig_dir,tag),w = 1000, h = 1100)
    layout(matrix(c(1,0,2,3),nrow= 2,ncol = 2),widths = c(400,500),heights = c(900,80))
    
    
    par(mar= c(0,4,8,6))
    image(pmin(pmax(t(lfp_emb[gnms_f,]),-2),2),col = shades,breaks = breaks,axes = F)
    mtext(text = substr(gnms_f,1,9),side = 4,at = seq(0,1,length.out = length(gnms_f)),las = 2)
    mtext(text = "Diff. expression per embryo \n KO vs wt for embryonic mesoderm",side = 3,at = 0.5,cex = 2)
    mtext(text = colnames(lfp_emb),side = 1,at = seq(0,1,length.out = ncol(lfp_emb)),las = 2)
    
    
    par(mar= c(0,4,8,6))
    image(pmin(pmax(t(lfp_ct[gnms_f,]),-2),2),col = shades,breaks = breaks,axes = F)
    mtext(text = substr(gnms_f,1,9),side = 4,at = seq(0,1,length.out = length(gnms_f)),las = 2)
    mtext(text = "Diff. expression KO vs wt \n per cell type",side = 3,at = 0.5,cex = 2)
    
    par(mar = c(4,4,0,6))
    image(as.matrix(c(1:ncol(lfp_ct))),col = colnames(lfp_ct),axes = F)
    
    
    dev.off()
    
  }
}



# generates the list of genes displayed in Figure 6E abd S7C
select_differentially_expressed_genes = function() {
  
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  mat_wt = scdb_mat("sing_emb_wt10")
  
  bad_genes = read.table(file = "data/external_data/sing_emb_wt10.bad_genes.txt",sep = "\t",stringsAsFactors = F)$x
  bad_genes = union(bad_genes,c("Igf2","AK145379;H19"))
  reg = 5e-5
  
  min_lfp = 0.7
  
  diff_lfp_meso_non_meso = 0.5
  
  mat_nm = "foxc_chim_wt10"
  mat = scdb_mat(mat_nm)
  
  # Caudal, paraxial and rostral mesoderm
  included_cts = c("#1a3f52","#408DA1","#8DB5CE")
  
  load(file = sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
  
  cls_f = colnames(mat@mat)[mat@cell_metadata[colnames(mat@mat),"cell_type"] == "KO"]
  cls_f = cls_f[cls_f %in% names(cmp_annot$query_cls_col)]
  cls_f = cls_f[cmp_annot$query_cls_col[cls_f] %in% included_cts]
  cls_wt = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% included_cts]
  
  umi_mat = mat@mat[,cls_f]
  egc_query = rowSums(umi_mat[rownames(mc_wt@e_gc),])/sum(umi_mat[rownames(mc_wt@e_gc),])
  egc_ref = rowSums(mat_wt@mat[rownames(mc_wt@e_gc),cls_wt])/sum(mat_wt@mat[rownames(mc_wt@e_gc),cls_wt])
  
  lfp = log2(egc_query + reg) - log2(egc_ref + reg)
  
  # filter genes that are differentially expressed also in ectoderm, endoderm and non-embryonic mesoderm
  bulk_lfp =  foxc_bulk_lfp_endo_ecto_exe_meso(reg = reg)
  
  # candidate_genes - log fold change larger than min_lfp
  gnms_f = names(lfp)
  gnms_f = setdiff(gnms_f,bad_genes)  
  f1 = abs(lfp[gnms_f]) > min_lfp 
  gnms_f = gnms_f[f1]
  f2 = abs(bulk_lfp[gnms_f] - lfp[gnms_f]) >  diff_lfp_meso_non_meso
  gnms_f = gnms_f[f2] 
  
  if(!dir.exists("data/chimera_tetraploid_analysis")) {
    dir.create("data/chimera_tetraploid_analysis")
  }
  if(0) {
    write.table(gnms_f,file = "data/chimera_tetraploid_analysis/fig_6e_differentially_expressed_genes.txt",sep = "\t")
  }
  return(gnms_f)
}

# the following function is only used in select_differentially_expressed_genes() above
# it calculates bulk differential expression between chimera Foxc DKO cells and WT atlas cells
# from endoderm, ectoderm and extraembryonic mesoderm cell types
foxc_bulk_lfp_endo_ecto_exe_meso = function(reg = 1e-5) {
  
  mat_nm = "foxc_chim_wt10"
  mat = scdb_mat(mat_nm)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  mat_wt = scdb_mat("sing_emb_wt10")
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  
  excluded_cts = c("#C594BF","#DFCDE4","#1a3f52","#408DA1","#8DB5CE","#45d1c5","#53f1fc","#FBBE92","#EF5A9D")
  
  load(file = sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
  
  cls_f = colnames(mat@mat)[mat@cell_metadata[colnames(mat@mat),"cell_type"] == "KO"]
  cls_f = cls_f[cls_f %in% names(cmp_annot$query_cls_col)]
  cls_f = cls_f[!(cmp_annot$query_cls_col[cls_f] %in% excluded_cts)]
  
  cols_f = table(cmp_annot$query_cls_col[cls_f])
  cols_f = names(cols_f)[cols_f > 20]
  print(col_to_ct[cols_f])
  cls_f = cls_f[cmp_annot$query_cls_col[cls_f] %in% cols_f]
  cls_wt = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% cols_f]
  
  umi_mat = mat@mat[rownames(mc_wt@e_gc),cls_f]
  
  egc_query = t(tgs_matrix_tapply(umi_mat,cmp_annot$query_cls_col[cls_f],sum))
  rownames(egc_query) = rownames(mc_wt@e_gc)
  egc_query = t(t(egc_query)/colSums(egc_query))
  egc_ref = t(tgs_matrix_tapply(mat_wt@mat[rownames(mc_wt@e_gc),cls_wt],mc_wt@colors[mc_wt@mc[cls_wt]],sum))
  rownames(egc_ref) = rownames(mc_wt@e_gc)
  egc_ref = t(t(egc_ref)/colSums(egc_ref))
  
  egc_query = rowMeans(egc_query)
  egc_ref = rowMeans(egc_ref)
  
  lfp  = log2(egc_query + reg) - log2(egc_ref + reg)
  
  return(lfp)
}


  
cmp_egc_per_chimera_embryo = function(mat_nm,gnms_f = NULL,included_colors = NULL,reg = 1e-5) {
    # calculates bulk expression per chim embryo for host and ko cells for the included_colors 
    # calculates corresponding average bulk expression for the best matched wild type embryo 
    # returns e_gc_ls bulk expression for host and KO cells
    
    #load(file = "data/paper_data/foxc12_chimera/nn_sc_color.Rda")
    
    mat=  scdb_mat(mat_nm)
    mc_wt = scdb_mc("sing_emb_wt10_recolored")
    
    if(is.null(included_colors)) {
      included_colors = mc_wt@color_key$color
    }
    
    df_chim = read.table(sprintf("data/chimera_tetraploid_analysis/%s/time_match/time_match_summary.txt",mat_nm),sep = "\t",h = T,stringsAsFactors = F)
    rownames(df_chim) = df_chim$embryo
    chim_embryos = df_chim$embryo
    
    # read single cell color information
    load(sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
    query_cls_col = cmp_annot$query_cls_col
    
    host_cls = names(query_cls_col)[(mat@cell_metadata[names(query_cls_col),"cell_type"] == "host") & 
                                      (mat@cell_metadata[names(query_cls_col),"embryo"] %in% chim_embryos)& 
                                      (query_cls_col %in% included_colors)] 
    query_cls = names(query_cls_col)[(mat@cell_metadata[names(query_cls_col),"cell_type"] %in% c("control","KO")) & 
                                       (mat@cell_metadata[names(query_cls_col),"embryo"] %in% chim_embryos) & 
                                       (query_cls_col %in% included_colors) ] 
    
    if(is.null(gnms_f)) {
      gnms_f = rownames(mat@mat)
    }
    
    egc_query = t(tgs_matrix_tapply(mat@mat[gnms_f,query_cls],mat@cell_metadata[query_cls,"embryo"],sum))
    rownames(egc_query) = gnms_f
    n_umis = colSums(mat@mat[,query_cls])
    n_umis = tapply(n_umis,mat@cell_metadata[query_cls,"embryo"],sum)
    
    egc_query = t(t(egc_query)/as.numeric(n_umis[colnames(egc_query)]))
    
    egc_host = t(tgs_matrix_tapply(mat@mat[gnms_f,host_cls],mat@cell_metadata[host_cls,"embryo"],sum))
    rownames(egc_host) = gnms_f
    n_umis = colSums(mat@mat[,host_cls])
    n_umis = tapply(n_umis,mat@cell_metadata[host_cls,"embryo"],sum)
    
    egc_host = t(t(egc_host)/as.numeric(n_umis[colnames(egc_host)]))
    
    legc_query = log2(egc_query + reg)
    legc_host = log2(egc_host + reg)
    
    return(list(egc_query = egc_query,egc_host = egc_host,legc_query = legc_query,legc_host = legc_host))
}
  
cmp_egc_wt_per_embryo = function(gnms_f = NULL,included_colors = NULL,reg = 1e-5) {
    
    mat = scdb_mat("sing_emb_wt10")
    mc_wt = scdb_mc("sing_emb_wt10_recolored")
    
    if(is.null(gnms_f)) {
      gnms_f = rownames(mat@mat)
    }
    
    if(is.null(included_colors)) {
      included_colors = mc@color_key$color
    }
    
    wt10_cls = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% included_colors]
    
    egc_wt = t(tgs_matrix_tapply(mat@mat[gnms_f,wt10_cls],mat@cell_metadata[wt10_cls,"transcriptional_rank"],sum))
    rownames(egc_wt) = gnms_f
    n_umis = colSums(mat@mat[,wt10_cls])
    n_umis = tapply(n_umis,mat@cell_metadata[wt10_cls,"transcriptional_rank"],sum)
    
    egc_wt = t(t(egc_wt)/as.numeric(n_umis[colnames(egc_wt)]))
    
    legc_wt = log2(egc_wt + reg)
    
    return(legc_wt)
}
  


if(0) {
  
  foxc_bulk_lfp = function(reg = 1e-5) {
    
    mat_nm = "foxc_chim_wt10"
    mat = scdb_mat(mat_nm)
    mc_wt = scdb_mc("sing_emb_wt10_recolored")
    mat_wt = scdb_mat("sing_emb_wt10")
    
    excluded_cts = c("#C594BF","#DFCDE4","#1a3f52","#408DA1","#8DB5CE","#45d1c5","#53f1fc","#FBBE92","#EF5A9D")
    
    load(file = sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
    
    cls_f = colnames(mat@mat)[mat@cell_metadata[colnames(mat@mat),"cell_type"] == "KO"]
    cls_f = cls_f[cls_f %in% names(cmp_annot$query_cls_col)]
    cls_f = cls_f[!(cmp_annot$query_cls_col[cls_f] %in% excluded_cts)]
    
    cols_f = table(cmp_annot$query_cls_col[cls_f])
    cols_f = names(cols_f)[cols_f > 20]
    print(cols_f)
    cls_f = cls_f[cmp_annot$query_cls_col[cls_f] %in% cols_f]
    
    umi_mat = mat@mat[rownames(mc_wt@e_gc),cls_f]
    
    egc_query = tgs_matrix_tapply(umi_mat,cmp_annot$query_cls_col[cls_f],sum)
    egc_query = t(t(egc_query)/colSums(egc_query))
    egc_re
    
    
    cmp = proj_on_wt_atlas(umi_mat = umi_mat,reg = reg)
    
    return(cmp)
  }
  
  
  proj_on_wt_atlas = function(umi_mat,ignore_genes = NULL,reg = 1e-5) {
    
    gset = scdb_gset("sing_emb_wt10")
    feat_genes = names(gset@gene_set)
    
    if(!is.null(ignore_genes)) {
      feat_genes = setdiff(feat_genes,ignore_genes)
    }
    
    mc_wt = scdb_mc("sing_emb_wt10_recolored")
    mat_wt = scdb_mat("sing_emb_wt10")
    
    n_umi_mc = t(tgs_matrix_tapply(mat_wt@mat[rownames(mc_wt@e_gc),names(mc_wt@mc)],mc_wt@mc,sum))
    rownames(n_umi_mc) = rownames(mc_wt@e_gc)
    
    sd_umi_mc = sqrt(n_umi_mc)
    sd_umi_mc = t(t(sd_umi_mc)/colSums(n_umi_mc))
    
    egc = t(t(n_umi_mc)/colSums(n_umi_mc))
    
    #egc = t(t(mc_wt@e_gc)/colSums(mc_wt@e_gc))
    
    legc = log2(egc + reg)
    
    
    
    best_cor = tgs_cor(as.matrix(umi_mat[feat_genes,]),legc[feat_genes,])
    
    best_ref = apply(best_cor,1,which.max)
    
    mat_n = t(t(umi_mat[rownames(mc_wt@e_gc),])/colSums(umi_mat[rownames(mc_wt@e_gc),]))
    
    egc_1 = rowMeans(mat_n)
    egc_2 = rowMeans(egc[,best_ref])
    
    sd_ref = rowMeans(sd_umi_mc[,best_ref])
    
    lfp = log2(egc_1 + reg) - log2(egc_2 + reg)
    return(list(egc_query = egc_1,egc_ref = egc_2, lfp = lfp,ref_mat = egc[,best_ref],query_mat = mat_n,best_ref_mc = best_ref,sd_ref = sd_ref))
  }
  
  
  

  
}
