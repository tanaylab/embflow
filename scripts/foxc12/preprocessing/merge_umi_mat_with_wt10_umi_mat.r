library("tidyverse")

merge_umi_mat_with_wt10 = function(mat_nm,new_mat_nm) {
  mat_nm1 = mat_nm
  mat_nm2 = "sing_emb_wt10"
  
  new_mat_id = new_mat_nm
  
  mat1 = scdb_mat(mat_nm1)
  mat2 = scdb_mat(mat_nm2)
  
  mat1@cell_metadata$cell_type = as.character(mat1@cell_metadata$cell_type)
  
  ignored_cls = union(mat1@ignore_cells,mat2@ignore_cells)
  ignored_genes = union(mat1@ignore_genes,mat2@ignore_genes)
  
  mat_ls = c(mat1,mat2)
  
  mat_all = rbind(cbind(mat1@mat,mat1@ignore_cmat),cbind(mat1@ignore_gmat,mat1@ignore_gcmat))
  md_all = mat1@cell_metadata
  
  
  i = 2
  mat = mat_ls[[i]]
  mat_tmp = rbind(cbind(mat@mat,mat@ignore_cmat),cbind(mat@ignore_gmat,mat@ignore_gcmat))
  mat_tmp = mat_tmp[rownames(mat_all),]
  mat_all = cbind(mat_all,mat_tmp)
  md_all = bind_rows(md_all,mat@cell_metadata)
  
  rownames(md_all) = md_all$cell
  md_all[colnames(mat2@mat),"cell_type"] = "wt10"
  
  
  
  mat_new = scm_new_matrix(mat = mat_all, stat_type =  "umi",cell_metadata = md_all)
  mat_new = scm_ignore_cells(scmat = mat_new,ig_cells = ignored_cls)
  mat_new = scm_ignore_genes(scmat = mat_new,ig_genes = ignored_genes)
  
  
  scdb_add_mat(id = new_mat_id,mat = mat_new)
  
}
