library("dplyr")

arsinh = function(x,loc = 0,w = 1) {
  
  return(log((x - loc)/w + sqrt(((x - loc)/w)^2 + 1)))
}


calc_a_b_from_x_y = function(x1,y1,x2,y2) {

  b = (y2 - y1)/(x2 - x1)
  a = (y1*x2 - x1*y2)/(x2 - x1) 
  return(c(a,b))
}

gating_of_foxc_control = function(fig_dir) {
  
  mat_nm = "foxc_control"
  mat = scdb_mat(mat_nm)
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  cls = mat@cell_metadata$cell[mat@cell_metadata$embryo != "empty"]
  
  cls_type = rep("empty",nrow(mat@cell_metadata))
  names(cls_type) = rownames(mat@cell_metadata)
  
  cls_type[cls] = "unclear"
  
  # parameters needed for the transformation
  loc_x = -600
  w_x = 400
  
  gfp_tr = arsinh(mat@cell_metadata[cls,"gfp_a"], loc = loc_x,w = w_x)
  ssc_w = mat@cell_metadata[cls,"ssc_w"]
  fsc_a = mat@cell_metadata[cls,"fsc_a"]
  xlim = c(min(gfp_tr),max(gfp_tr))
  ylim = c(min(fsc_a),max(fsc_a))
  
  col_fsc_ssc = densCols(x = gfp_tr,y = fsc_a)
  
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "gfp",col = col_fsc_ssc,
       xlim = xlim,ylim = ylim)
  
  names(gfp_tr) = cls
  names(ssc_w) = cls
  names(fsc_a) = cls
  
  pdf(sprintf("%s/gfp_vs_fsc_all_cls.pdf",fig_dir),useDingbats = F)
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "GFP",main = "FoxC and control chimera and tetraploid cells",
       xlim = xlim,ylim = ylim)
  dev.off()
  
  # Gating of control chimera
  
  loc_x = -300
  w_x = 400
  
  a_l = -180000
  b_l = 150000
  a_h = -290000
  b_h = 150000
  
  
  cls_f = cls[mat@cell_metadata[cls,"Experiment_type"] %in% c("Control chimera")]
  
  ssc_w = mat@cell_metadata[cls_f,"ssc_w"]
  fsc_a = mat@cell_metadata[cls_f,"fsc_a"]
  gfp_tr = arsinh(mat@cell_metadata[cls_f,"gfp_a"], loc = loc_x,w = w_x)
  names(gfp_tr) = cls_f
  names(fsc_a) = cls_f
  xlim = c(min(gfp_tr),max(gfp_tr))
  #ylim = c(min(fsc_a),max(fsc_a))
  
  col_dens = densCols(gfp_tr,fsc_a)
  pdf(sprintf("%s/control_chimera_gfp_vs_fsc.pdf",fig_dir),useDingbats = F)
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "gfp",col= col_dens,main = "Control chimera cells",
       xlim = xlim,ylim = ylim)
  
  abline(a = a_l,b = b_l,lty = "dashed")
  abline(a = a_h,b = b_h,lty = "dashed")
  dev.off()
  
  
  f_control = fsc_a < a_h + gfp_tr*b_h
  f_host = fsc_a > a_l + gfp_tr*b_l
  
  col_cls= rep("gray",length(cls_f))
  col_cls[f_control] = "cornflowerblue"
  col_cls[f_host] = "black"
  
  pdf(sprintf("%s/control_chimera_gfp_vs_fsc_gated.pdf",fig_dir),useDingbats = F)
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "GFP",main = "Control chimera cells",col = col_cls,
       xlim = xlim,ylim = ylim)
  
  abline(a = a_l,b = b_l,lty = "dashed")
  abline(a = a_h,b = b_h,lty = "dashed")
  
  legend(x = "topleft",legend = c("unclear","control","Host"),col = c("gray","cornflowerblue","black"),pch = 19)
  dev.off()
  
  
  cls_type[cls_f[f_control]] = "control"
  cls_type[cls_f[f_host]] = "host"
  
  embryos = unique(mat@cell_metadata[cls_f,"embryo"])

  
  # Gating of control tetraploids
  
  loc_x = -300
  w_x = 400
  
  a_l = -180000
  b_l = 150000
  a_h = -290000
  b_h = 150000
  
  
  cls_f = cls[mat@cell_metadata[cls,"Experiment_type"] %in% c("Control tetraploid")]
  
  ssc_w = mat@cell_metadata[cls_f,"ssc_w"]
  fsc_a = mat@cell_metadata[cls_f,"fsc_a"]
  gfp_tr = arsinh(mat@cell_metadata[cls_f,"gfp_a"], loc = loc_x,w = w_x)
  names(gfp_tr) = cls_f
  names(fsc_a) = cls_f
  xlim = c(min(gfp_tr),max(gfp_tr))
  #ylim = c(min(fsc_a),max(fsc_a))
  
  col_dens = densCols(gfp_tr,fsc_a)
  pdf(sprintf("%s/control_tetra_gfp_vs_fsc.pdf",fig_dir),useDingbats = F)
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "gfp",col= col_dens,main = "Control tetraploid cells",
       xlim = xlim,ylim = ylim)
  
  abline(a = a_l,b = b_l,lty = "dashed")
  abline(a = a_h,b = b_h,lty = "dashed")
  dev.off()
  
  
  f_control = fsc_a < a_h + gfp_tr*b_h
  f_host = fsc_a > a_l + gfp_tr*b_l
  
  col_cls= rep("gray",length(cls_f))
  col_cls[f_control] = "cornflowerblue"
  col_cls[f_host] = "black"
  
  pdf(sprintf("%s/control_tetra_gfp_vs_fsc_gated.pdf",fig_dir),useDingbats = F)
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "GFP",main = "Control tetraploid cells",col = col_cls,
       xlim = xlim,ylim = ylim)
  
  abline(a = a_l,b = b_l,lty = "dashed")
  abline(a = a_h,b = b_h,lty = "dashed")
  
  legend(x = "topleft",legend = c("unclear","control","Host"),col = c("gray","cornflowerblue","black"),pch = 19)
  dev.off()
  
  
  cls_type[cls_f[f_control]] = "control"
  cls_type[cls_f[f_host]] = "host"
  
  embryos = unique(mat@cell_metadata[cls_f,"embryo"])
  

  # Gating of Foxc chimera
  
  loc_x = 0
  w_x = 50
  
  a_l = -60000
  b_l = 60000
  a_h = -110000
  b_h = 60000
  
  
  cls_f = cls[mat@cell_metadata[cls,"Experiment_type"] %in% c("Foxc chimera")]
  
  fsc_a = mat@cell_metadata[cls_f,"fsc_a"]
  ssc_w = mat@cell_metadata[cls_f,"ssc_w"]
  gfp_tr = arsinh(mat@cell_metadata[cls_f,"gfp_a"], loc = loc_x,w = w_x)
  names(gfp_tr) = cls_f
  names(fsc_a) = cls_f
  xlim = c(min(gfp_tr),max(gfp_tr))
  #ylim = c(min(fsc_a),max(fsc_a))
  
  col_dens = densCols(gfp_tr,fsc_a)
  pdf(sprintf("%s/foxc_chimera_gfp_vs_fsc.pdf",fig_dir),useDingbats = F)
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "GFP",col= col_dens,main = "FoxC chimera cells",
       xlim = xlim,ylim = ylim)
  
  abline(a = a_l,b = b_l,lty = "dashed")
  abline(a = a_h,b = b_h,lty = "dashed")
  dev.off()
  
  
  f_ko = fsc_a < a_h + gfp_tr*b_h
  f_host = fsc_a > a_l + gfp_tr*b_l
  
  col_cls= rep("gray",length(cls_f))
  col_cls[f_ko] = "#83c26d"
  col_cls[f_host] = "black"
  
  pdf(sprintf("%s/foxc_chimera_gfp_vs_fsc_gated.pdf",fig_dir),useDingbats = F)
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "GFP",main = "FoxC chimera cells",col = col_cls,
       xlim = xlim,ylim = ylim)
  
  abline(a = a_l,b = b_l,lty = "dashed")
  abline(a = a_h,b = b_h,lty = "dashed")
  
  legend(x = "topleft",legend = c("unclear","KO","Host"),col = c("gray","#83c26d","black"),pch = 19)
  dev.off()
  
  
  cls_type[cls_f[f_ko]] = "KO"
  cls_type[cls_f[f_host]] = "host"
  
  embryos = unique(mat@cell_metadata[cls_f,"embryo"])
  

  
  # Gating of Foxc tetraploid
  
  loc_x = -700
  w_x = 400
  
  a_l = -200000
  b_l = 100000
  a_h = -260000
  b_h = 100000
  
  cls_f = cls[mat@cell_metadata[cls,"Experiment_type"] %in% c("Foxc tetraploid")]
  
  fsc_a = mat@cell_metadata[cls_f,"fsc_a"]
  ssc_w = mat@cell_metadata[cls_f,"ssc_w"]
  gfp_tr = arsinh(mat@cell_metadata[cls_f,"gfp_a"], loc = loc_x,w = w_x)
  names(gfp_tr) = cls_f
  names(fsc_a) = cls_f
  xlim = c(min(gfp_tr),max(gfp_tr))
  #ylim = c(min(fsc_a),max(fsc_a))
  
  col_dens = densCols(gfp_tr,fsc_a)
  pdf(sprintf("%s/foxc_tetraploid_gfp_vs_fsc.pdf",fig_dir),useDingbats = F)
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "gfp",col= col_dens,main = "FoxC tetraploid cells",
       xlim = xlim,ylim = ylim)
  
  abline(a = a_l,b = b_l,lty = "dashed")
  abline(a = a_h,b = b_h,lty = "dashed")
  dev.off()
  
  
  f_ko = fsc_a < a_h + gfp_tr*b_h
  f_host = fsc_a > a_l + gfp_tr*b_l
  
  col_cls= rep("gray",length(cls_f))
  col_cls[f_ko] = "#83c26d"
  col_cls[f_host] = "black"
  
  pdf(sprintf("%s/foxc_tetraploid_gfp_vs_fsc_gated.pdf",fig_dir),useDingbats = F)
  plot(gfp_tr,fsc_a,pch = 19,cex = 0.4,ylab = "FSC-A",xlab = "GFP",main = "FoxC tetraploid cells",col = col_cls,
       xlim = xlim,ylim = ylim)
  
  abline(a = a_l,b = b_l,lty = "dashed")
  abline(a = a_h,b = b_h,lty = "dashed")
  
  legend(x = "topleft",legend = c("unclear","KO","Host"),col = c("gray","#83c26d","black"),pch = 19)
  dev.off()
  
  
  cls_type[cls_f[f_ko]] = "KO"
  cls_type[cls_f[f_host]] = "host"
  
  embryos = unique(mat@cell_metadata[cls_f,"embryo"])
  
  
  df_gating = data.frame(cell = rownames(mat@cell_metadata),cell_type = cls_type,stringsAsFactors = F)
  
  if(0) {
    
    df_gating = data.frame(cell = rownames(mat@cell_metadata),cell_type = cls_type,stringsAsFactors = F)
    mat = scdb_mat(mat_nm)
    md = mat@cell_metadata
    md = left_join(md,df_gating,by = "cell")
    rownames(md) = md$cell
    mat@cell_metadata = md
    scdb_add_mat(id = mat_nm,mat = mat)
  }
  
  return(df_gating)
}


