
transfer_color_chimera_tetraploid = function(mat_nm,tag = "query") {
  
  cgraph_id = mat_nm
  mat = scdb_mat(mat_nm)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  
  host_cls = colnames(mat@mat)[( mat@cell_metadata[colnames(mat@mat),"cell_type"] == "host" )]
  injected_cls = colnames(mat@mat)[( mat@cell_metadata[colnames(mat@mat),"cell_type"] %in% c("KO","control") )]
  wt10_cls = names(mc_wt@mc)
  
  
  query_cls = c(host_cls,injected_cls)
  query_cls_md = c(rep("host",length(host_cls)),rep(tag,length(injected_cls)))
  names(query_cls_md) = query_cls
  
  ref_cls_color = mc_wt@colors[mc_wt@mc[wt10_cls]]
  names(ref_cls_color)  = wt10_cls
  
  cmp_annot = color_query_by_reference(cgraph_id = cgraph_id,query_cls_md = query_cls_md,ref_cls_color = ref_cls_color)
  
  data_dir = sprintf("data/chimera_tetraploid_analysis/%s",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  data_dir = sprintf("data/chimera_tetraploid_analysis/%s/color_annotation",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  print(data_dir)
  save(cmp_annot,file = sprintf("%s/cmp_annot.Rda",data_dir))
  
}




color_query_by_reference = function(cgraph_id,
                                    query_cls_md,
                                    ref_cls_color,
                                    max_knn_cgraph = 50,
                                    knn_color_space = 1,
                                    color_query_cls_with_few_ref_neighbors = T,
                                    threshold_few_ref_neighbors = 0.5) {
  
  cgraph = scdb_cgraph(cgraph_id)
  
  query_cls_f = intersect(names(query_cls_md),cgraph@nodes)
  
  ref_cls = names(ref_cls_color)
  ref_cls_f = intersect(ref_cls,cgraph@nodes)
  all_cls = c(query_cls_f,ref_cls_f)
  
  cls1 = levels(cgraph@edges$mc1)
  cls2 = levels(cgraph@edges$mc2)
  
  # add color and type to 
  cl_to_type = c(query_cls_md[query_cls_f],rep("ref",length(ref_cls_f)))
  names(cl_to_type) = all_cls
  
  type_cls1 = cl_to_type[cls1]
  names(type_cls1) = cls1
  type_cls2 = cl_to_type[cls2]
  names(type_cls2) = cls2
  
  color_cls1 = ref_cls_color[cls1]
  names(color_cls1) = cls1
  color_cls2 = ref_cls_color[cls2]
  names(color_cls2) = cls2
  
  cgraph@edges$type1 = type_cls1[cgraph@edges$mc1]
  cgraph@edges$type2 = type_cls2[cgraph@edges$mc2]  
  cgraph@edges$color1 = color_cls1[cgraph@edges$mc1]
  cgraph@edges$color2 = color_cls2[cgraph@edges$mc2]
  
  
  # report on cells which have few reference neighbors
  cl_vs_type = table(cgraph@edges$mc1,cgraph@edges$type2)
  cl_vs_type = cl_vs_type[all_cls,]
  cls_with_few_neighbors = rownames(cl_vs_type)[rowSums(cl_vs_type) <10]
  cl_vs_type_n = cl_vs_type[rowSums(cl_vs_type) > 0,]
  cl_vs_type_n = cl_vs_type_n/rowSums(cl_vs_type_n)
  
  cls_with_low_fraction_of_ref_cell_neighbors = rownames(cl_vs_type_n)[1 - cl_vs_type_n[,"ref"] > threshold_few_ref_neighbors]
  
  if(!color_query_cls_with_few_ref_neighbors) {
    query_cls_f = setdiff(query_cls_f,cls_with_low_fraction_of_ref_cell_neighbors)
  }
  
  # color annotation of query_cls - restrict to max_knn_cgraph neighbors in the cgraph
  
  f = cgraph@edges$w > (1 - max_knn_cgraph/100)
  graph_f = cgraph@edges[f,]
  
  
  # calculate color contingency table for query vs ref and ref vs ref
  cl_vs_color = table(graph_f$mc1,graph_f$color2)
  cl_vs_color_ref = cl_vs_color[ref_cls_f,]
  cl_vs_color_query = cl_vs_color[query_cls_f,]
  cl_vs_color_ref = cl_vs_color_ref[rowSums(cl_vs_color_ref) > 0,]
  cl_vs_color_query = cl_vs_color_query[rowSums(cl_vs_color_query) > 0,]
  
  ref_ref_col_cor = tgs_cor_knn(x = t(cl_vs_color_ref),y = t(cl_vs_color_ref),knn = (knn_color_space + 1))
  query_ref_col_cor = tgs_cor_knn(x = t(cl_vs_color_query),y = t(cl_vs_color_ref),knn = knn_color_space)
  
  
  sample_nn = sample(2:(knn_color_space + 1),size = nrow(cl_vs_color_ref),replace = T)
  sample_ind = c(0:(nrow(cl_vs_color_ref)-1))*(knn_color_space + 1) + sample_nn
  
  original_cl = ref_ref_col_cor$col1[sample_ind]
  ref_sampled_color = ref_cls_color[as.character(ref_ref_col_cor$col2[sample_ind])]
  ref_original_color = ref_cls_color[as.character(original_cl)]
  names(ref_original_color) = original_cl
  names(ref_sampled_color) = original_cl
  
  
  # sample color for query cls
  sample_nn = sample(1:knn_color_space,size = nrow(cl_vs_color_query),replace = T)
  sample_ind = c(0:(nrow(cl_vs_color_query)-1))*knn_color_space + sample_nn
  
  query_cls_col = ref_cls_color[as.character(query_ref_col_cor$col2[sample_ind])]
  names(query_cls_col) = as.character(query_ref_col_cor$col1[sample_ind])
  
  cmp_col_annot = list(query_cls_col = query_cls_col,
                       cls_with_low_fraction_of_ref_cell_neighbors = cls_with_low_fraction_of_ref_cell_neighbors,
                       ref_cls_orig_color = ref_original_color,
                       ref_cls_sampled_color = ref_sampled_color,
                       cls_with_few_neighbors = cls_with_few_neighbors,
                       max_knn_cgraph = max_knn_cgraph,
                       knn_color_space = knn_color_space,
                       query_cls_col_dist = cl_vs_color_query,
                       ref_cls_col_dist = cl_vs_color_ref,
                       cl_type_dist = cl_vs_type)
  
  return(cmp_col_annot)
}







