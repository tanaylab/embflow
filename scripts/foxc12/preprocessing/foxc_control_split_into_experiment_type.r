source("scripts/foxc12/preprocessing/merge_umi_mat_with_wt10_umi_mat.r")

if(0) {
  
  ignore_embryos = c("6D_m2e3","6D_m3e1","6D_m3e7")
  
  mat_nm = "foxc_control"
  mat = scdb_mat(mat_nm)
  
  
  df_nm = data.frame(type = c("Control chimera","Control tetraploid","Foxc chimera","Foxc tetraploid"),
                     mat_nm = c("control_chim","control_tetra","foxc_chim","foxc_tetra"),stringsAsFactors = F)
  
  for (i in 1:4) {
    
    exp_nm = df_nm[i,"type"]
    new_mat_nm = df_nm[i,"mat_nm"]
    
    cls = colnames(mat@mat)[mat@cell_metadata[colnames(mat@mat),"Experiment_type"] == exp_nm]
    cls = cls[!(mat@cell_metadata[cls,"embryo"] %in% ignore_embryos)]
      
    mat_new = scm_ignore_cells(scmat = mat,ig_cells = cls,reverse = T)
    
    scdb_add_mat(id = new_mat_nm,mat = mat_new)
    
    merge_umi_mat_with_wt10(mat_nm = new_mat_nm,new_mat_nm = paste0(new_mat_nm,"_wt10"))
  }
  
  
}