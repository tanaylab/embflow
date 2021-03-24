library("metacell")
scdb_init("scrna_db/",force_reinit = T)
tgconfig::override_params(config_file = "config/sing_emb.yaml",package = "metacell")
source("scripts/foxc12/transfer_cell_type_annotation.r")
source("scripts/foxc12/transfer_time_annotation.r")
source("scripts/foxc12/foxc_chim_timing.r")
source("scripts/foxc12/control_chim_timing.r")
source("scripts/foxc12/definition_cell_types_endoderm_ectoderm_embryonic_mesoderm.r")

foxc_chimera_generate_time_and_cell_type_annotation = function() {
  
  if(!dir.exists("data/chimera_tetraploid_analysis")) {
    dir.create("data/chimera_tetraploid_analysis")
  }
  # generate single-cell graph
  gen_cgraph(mat_nm = "foxc_chim_wt10")
  message("generated cgraph")
  transfer_color_chimera_tetraploid(mat_nm = "foxc_chim_wt10",tag = "KO")
  message("transfered cell type annotation")
  foxc_chim_timing()
  message("tranfered time annotation")
  
  # generate time distributions per embryo based on endo- and ectoderm cells
  # endo_ct_colors() and ectoderm_ct_colors just return the colors of ecto- and endodermal cell types
  endo_colors = endo_ct_colors()
  ecto_colors = ectoderm_ct_colors()
  included_colors = c(endo_colors,ecto_colors)
  foxc_chim_timing(tag = "endo_ecto",included_colors = included_colors)
  
  # Embryonic mesoderm
  included_colors = embryonic_meso_ct_colors()
  foxc_chim_timing(tag = "emb_meso",included_colors = included_colors)
  
}

control_chimera_generate_time_and_cell_type_annotation = function() {
  
  if(!dir.exists("data/chimera_tetraploid_analysis")) {
    dir.create("data/chimera_tetraploid_analysis")
  }
  # generate single-cell graph
  gen_cgraph(mat_nm = "control_chim_wt10")
  message("generated cgraph")
  transfer_color_chimera_tetraploid(mat_nm = "control_chim_wt10",tag = "control")
  message("transfered cell type annotation")
  control_chim_timing()
  message("tranfered time annotation")
  
  # generate time distributions per embryo based on endo- and ectoderm cells
  # endo_ct_colors() and ectoderm_ct_colors just return the colors of ecto- and endodermal cell types
  endo_colors = endo_ct_colors()
  ecto_colors = ectoderm_ct_colors()
  included_colors = c(endo_colors,ecto_colors)
  control_chim_timing(tag = "endo_ecto",included_colors = included_colors)
  
  # Embryonic mesoderm
  included_colors = embryonic_meso_ct_colors()
  control_chim_timing(tag = "emb_meso",included_colors = included_colors)
  
}

foxc_tetraploid_generate_time_and_cell_type_annotation = function() {
  
  if(!dir.exists("data/chimera_tetraploid_analysis")) {
    dir.create("data/chimera_tetraploid_analysis")
  }
  # generate single-cell graph
  gen_cgraph(mat_nm = "foxc_tetra_wt10")
  message("generated cgraph")
  transfer_color_chimera_tetraploid(mat_nm = "foxc_tetra_wt10",tag = "KO")
  
  foxc_tetra_timing()
  
  
}

control_tetraploid_generate_time_and_cell_type_annotation = function() {
  
  if(!dir.exists("data/chimera_tetraploid_analysis")) {
    dir.create("data/chimera_tetraploid_analysis")
  }
  # generate single-cell graph
  gen_cgraph(mat_nm = "control_tetra_wt10")
  message("generated cgraph")
  transfer_color_chimera_tetraploid(mat_nm = "control_tetra_wt10",tag = "control")
  
  control_tetra_timing()
  
  
}


gen_cgraph = function(mat_nm,T_vm = 0.1,Knn = 100) {
  
  bad_genes = read.table(file = "data/external_data/sing_emb_wt10.bad_genes.txt",sep = "\t",stringsAsFactors = F)$x
  bad_genes = c(bad_genes,c("Igf2","AK145379;H19","Polg","Slc25a4","Peg10","Igf2as","AK086477;Sh3glb1"))
  
  mcell_add_gene_stat(mat_nm, mat_nm, force=T)
  
  mcell_gset_filter_varmean(mat_nm, mat_nm, T_vm=T_vm, force_new=T)
  mcell_gset_filter_cov(mat_nm, mat_nm, T_tot=50, T_top3=3)
  
  gset = scdb_gset(mat_nm)
  nms = names(gset@gene_set)
  #bad gene that will be removed from list of genes that helps to mark metacell 
  bad_g = c(grep("^Rpl",nms,v=T),grep("^Gm",nms,v=T),grep("Rps",nms,v=T))
  
  bad_g = c(bad_g, bad_genes)
  
  gset_f = gset_new_restrict_nms(gset=gset, bad_g, inverse=T, "feat filt")
  scdb_add_gset(mat_nm, gset_f)
  
  mcell_add_cgraph_from_mat_bknn(mat_id=mat_nm, 
                                 gset_id = mat_nm, 
                                 graph_id=mat_nm,
                                 K=Knn,
                                 dsamp=T)
  
}


foxc_tetra_timing = function() {
  
  mat_nm = "foxc_tetra_wt10"
  mc_id = "foxc_tetra_wt10_recolored"
  cgraph_id =  mat_nm
  mat = scdb_mat(mat_nm)
  mc = scdb_mc(mc_id)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  
  excluded_colors = c("#F6BFCB","#7F6874")
  
  load(sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
  
  f = mat@cell_metadata[colnames(mat@mat),"cell_type"] %in% c("KO")
  tetra_cls = colnames(mat@mat)[f]
  tetra_cls = intersect(tetra_cls,names(cmp_annot$query_cls_col))
  
  excluded_cls = tetra_cls[cmp_annot$query_cls_col[tetra_cls] %in% excluded_colors]
  tetra_cls = setdiff(tetra_cls,excluded_cls)
  
  ko_embryos = unique(mat@cell_metadata[tetra_cls,"embryo"])
  
  wt10_cls = intersect(names(mc_wt@mc)[ !(mc_wt@colors[mc_wt@mc] %in% excluded_colors) ],colnames(mat@mat))
  atlas_time = mat@cell_metadata[wt10_cls,"transcriptional_rank"]
  names(atlas_time) = wt10_cls
  
  data_dir = sprintf("data/chimera_tetraploid_analysis/%s",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  
  data_dir = sprintf("data/chimera_tetraploid_analysis/%s/time_match",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  
  
  
  # first get atlas time distribution
  
  atlas_time_dist = get_atlas_time_dist(atlas_time = atlas_time,graph_id = cgraph_id)
  
  # timing using KO cells
  
  ko_cls = tetra_cls[( mat@cell_metadata[tetra_cls,"cell_type"] == "KO" ) & ( mat@cell_metadata[tetra_cls,"embryo"] %in% ko_embryos )]
  query_cls_md = mat@cell_metadata[ko_cls,"embryo"]
  names(query_cls_md) = ko_cls
  
  query_time_dist_ko = get_query_time_dist(query_cls_md = query_cls_md,atlas_time = atlas_time,graph_id = cgraph_id)
  
  
  time_dist_ko = list(atlas_time_dist = atlas_time_dist$atlas_time_dist,
                      ko = query_time_dist_ko$query_time_dist)
  
  save(time_dist_ko,file = sprintf("%s/time_dist_ko.Rda",data_dir))
  
  chim_emb_summary = as.data.frame.matrix(table(mat@cell_metadata[tetra_cls,"embryo"],mat@cell_metadata[tetra_cls,"cell_type"]))
  chim_emb_summary$embryo = rownames(chim_emb_summary)
  chim_emb_summary$best_rank_ko = NA
  chim_emb_summary[rownames(time_dist_ko$ko),"best_rank_ko"] = time_dist_best_match(atlas_time_dist = time_dist_ko$atlas_time_dist,
                                                                                    query_time_dist = time_dist_ko$ko)
  
  write.table(chim_emb_summary,file = sprintf("%s/time_match_summary.txt",data_dir),sep ="\t",row.names = F)
}


control_tetra_timing = function() {
  
  mat_nm = "control_tetra_wt10"
  mc_id = "control_tetra_wt10_recolored"
  cgraph_id =  mat_nm
  mat = scdb_mat(mat_nm)
  mc = scdb_mc(mc_id)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  
  excluded_colors = c("#F6BFCB","#7F6874")
  
  load(sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
  
  f = mat@cell_metadata[colnames(mat@mat),"cell_type"] %in% c("control")
  tetra_cls = colnames(mat@mat)[f]
  tetra_cls = intersect(tetra_cls,names(cmp_annot$query_cls_col))
  
  excluded_cls = tetra_cls[cmp_annot$query_cls_col[tetra_cls] %in% excluded_colors]
  tetra_cls = setdiff(tetra_cls,excluded_cls)
  
  control_embryos = unique(mat@cell_metadata[tetra_cls,"embryo"])
  
  wt10_cls = intersect(names(mc_wt@mc)[ !(mc_wt@colors[mc_wt@mc] %in% excluded_colors) ],colnames(mat@mat))
  atlas_time = mat@cell_metadata[wt10_cls,"transcriptional_rank"]
  names(atlas_time) = wt10_cls
  
  data_dir = sprintf("data/chimera_tetraploid_analysis/%s",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  
  data_dir = sprintf("data/chimera_tetraploid_analysis/%s/time_match",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  
  
  
  # first get atlas time distribution
  
  atlas_time_dist = get_atlas_time_dist(atlas_time = atlas_time,graph_id = cgraph_id)
  
  # timing using control cells
  
  control_cls = tetra_cls[( mat@cell_metadata[tetra_cls,"cell_type"] == "control" ) & ( mat@cell_metadata[tetra_cls,"embryo"] %in% control_embryos )]
  query_cls_md = mat@cell_metadata[control_cls,"embryo"]
  names(query_cls_md) = control_cls
  
  query_time_dist_control = get_query_time_dist(query_cls_md = query_cls_md,atlas_time = atlas_time,graph_id = cgraph_id)
  
  
  time_dist_control = list(atlas_time_dist = atlas_time_dist$atlas_time_dist,
                           control = query_time_dist_control$query_time_dist)
  
  save(time_dist_control,file = sprintf("%s/time_dist_control.Rda",data_dir))
  
  chim_emb_summary = as.data.frame.matrix(table(mat@cell_metadata[tetra_cls,"embryo"],mat@cell_metadata[tetra_cls,"cell_type"]))
  chim_emb_summary$embryo = rownames(chim_emb_summary)
  chim_emb_summary$best_rank_control = NA
  chim_emb_summary[rownames(time_dist_control$control),"best_rank_control"] = time_dist_best_match(atlas_time_dist = time_dist_control$atlas_time_dist,
                                                                                                   query_time_dist = time_dist_control$control)
  
  write.table(chim_emb_summary,file = sprintf("%s/time_match_summary.txt",data_dir),sep ="\t",row.names = F)
}



