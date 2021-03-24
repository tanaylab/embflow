# removes extraembryonic ectoderm and parietal endoderm cells

if(0) {
  mat_nm = "foxc_control"
  
  mc_wt = scdb_mc("sing_emb_wt9_bs500f")
  
  mat_query = scdb_mat(mat_nm)
  scdb_add_mat(paste0(mat_nm,"_w_exe_ecto"),mat_query)
  gset = scdb_gset("sing_emb_wt9")
  feat_genes = names(gset@gene_set)
  
  egc_type = t(tgs_matrix_tapply(mc_wt@e_gc[feat_genes,],mc_wt@colors,mean))
  egc_type = egc_type[,colnames(egc_type) != "gray"]
  rownames(egc_type) = feat_genes
  
  legc = log2(egc_type + 1e-5)
  
  query_ref_cor = tgs_cor(as.matrix(mat_query@mat[feat_genes,]),egc_type)
  
  best_ref_type = colnames(query_ref_cor)[apply(query_ref_cor,1,which.max)]
  
  f_exe_ecto_parietal_endo = best_ref_type %in% c("#1A1A1A","#989898")
  
  cls_f = rownames(query_ref_cor)[f_exe_ecto_parietal_endo]
  
  mat_new = scm_ignore_cells(scmat = mat_query,ig_cells = union(cls_f,mat_query@ignore_cells))
  
  scdb_add_mat(id = mat_nm,mat = mat_new)
}
