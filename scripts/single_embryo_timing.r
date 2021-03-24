library("metacell")
library("tidyverse")
scdb_init("scrna_db/",force_reinit = T)
# calculate embryo time

gen_single_embryo_timing = function() {
  
  # external ranking of embryos containing morphology ranks and size measurements
  emb_age = read.table(file = "data/external_data/single_embryo_morphology_size_information.txt",sep = "\t",h = T)
  
  # generating intrinsic transcriptional ranking of embryos,"intrinsic ranking" if Figure 1
  intrinsic_ranking = wt10_gen_intrinsinc_ranking()
  
  # similarity matrix displayed in Figure 1B
  ee_similarity_mat = intrinsic_ranking$similarity_mat
  
  emb_age = left_join(emb_age,intrinsic_ranking$embryo_final_order,by = "embryo")
  
  
  # translating transcriptional ranks into embryonic developmental times 6.x to 8.x
  developmental_time = calculate_developmental_time_from_size_measurement(emb_age_time = emb_age) 
  
  emb_age = left_join(emb_age,developmental_time,by = "transcriptional_rank")

  # projecting embryos on reference gastrulation atlas Pijuan-Sala et al Nature (2019)
  ref_age = calculate_ref_gastru_atlas_age()
  
  emb_age = left_join(emb_age,ref_age,by = "embryo")
  
  return(list(single_embryo_time = emb_age,transcriptional_similarity_matrix = ee_similarity_mat))
}


calculate_developmental_time_from_size_measurement = function(emb_age_time) {
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")

  emb_age_time = emb_age_time[,c("embryo","transcriptional_rank","morphology_rank","area","mouse_type","age_group")]
  #emb_age_time = unique(mat@cell_metadata[names(mc@mc),c("embryo","transcriptional_rank","morphology_rank","area","mouse_type","age_group")])
  
  age_group_median_rank = tapply(emb_age_time$transcriptional_rank,emb_age_time$age_group,median)
  min_rank = round(age_group_median_rank[1])
  max_rank = round(age_group_median_rank[13])
  
  min_age = 6.5
  max_age = 8.1
  
  f = (emb_age_time$mouse_type != "ICR") & !(is.na(emb_age_time$area))
  emb_age_time = emb_age_time[f,]
  
  f2 = emb_age_time$morphology_rank > 100 & log2(emb_age_time$area) < 16.3
  
  fit_rank_vs_area = smooth.spline(x = emb_age_time[!f2,"transcriptional_rank"],log2(emb_age_time[!f2,"area"]),spar = 0.9)

  
  min_x = fit_rank_vs_area$x[1]
  max_x = fit_rank_vs_area$x[length(fit_rank_vs_area$x)]
  
  extrapol_fit_x = fit_rank_vs_area$x
  extrapol_fit_y = fit_rank_vs_area$y
  if(min_x > 1) {
    
    new_fit_y = (extrapol_fit_y[10] - extrapol_fit_y[1])/(extrapol_fit_x[10] - extrapol_fit_x[1])*(c(1:(min_x-1)) - min_x) + extrapol_fit_y[1]
    extrapol_fit_x = c(c(1:(min_x-1)),extrapol_fit_x)
    extrapol_fit_y = c(new_fit_y,extrapol_fit_y)
  }
  
  if(max_x < 153) {
    
    new_fit_y = (extrapol_fit_y[length(extrapol_fit_y)] - extrapol_fit_y[(length(extrapol_fit_y) - 1)])/(extrapol_fit_x[length(extrapol_fit_x)] - extrapol_fit_x[(length(extrapol_fit_x) - 1)])*(c(1:(153-max_x))) + extrapol_fit_y[length(extrapol_fit_y)]
    extrapol_fit_x = c(extrapol_fit_x,c((max_x + 1):153))
    extrapol_fit_y = c(extrapol_fit_y,new_fit_y)
  }
  
  area_interp = rep(0,153)
  
  for(n in 1:153) {
    if (n %in% extrapol_fit_x) {
      n_ind = which(extrapol_fit_x == n)
      area_interp[n] = extrapol_fit_y[n_ind]
    } else {
      
      n_l_ind = max(which(extrapol_fit_x < n))
      n_h_ind = min(which(extrapol_fit_x > n))
      
      n_l = extrapol_fit_x[n_l_ind]
      n_h = extrapol_fit_x[n_h_ind]
      
      t_interp = (n - n_l)/(n_h - n_l)*extrapol_fit_y[n_h_ind] + (n_h - n)/(n_h - n_l)*extrapol_fit_y[n_l_ind]
      
      area_interp[n] = t_interp
    }
  }
  
  
  infered_time = min_age + (area_interp - area_interp[min_rank])/(area_interp[max_rank] - area_interp[min_rank])*(max_age - min_age)
  
  embryonic_time = data.frame(transcriptional_rank = c(1:153),developmental_time = infered_time)
  
  return(embryonic_time)
}


calculate_ref_gastru_atlas_age = function() {
  
  atlas = mcell_gen_atlas(mat_id = "emb_gotg", 
                          mc_id = "emb_gotg_bs500f", 
                          gset_id  = "emb_gotg", 
                          mc2d_id= "emb_gotg_bs500f")
  
  mat_id = "sing_emb_wt10"
  mat_ref_id = "emb_gotg"
  mc_id = "sing_emb_wt10_recolored"
  
  # the following cell types are ignored because they show strong embryo-to-embryo fluctuations in there frequency
  list_of_ignored_colors = c("#1A1A1A","#989898","#7F6874","#c9a997","#C72228")
  
  ref_mc = scdb_mc(atlas$mc_id)
  query_mc = scdb_mc(mc_id)
  mat = scdb_mat(mat_id)
  mat_ref = scdb_mat(mat_ref_id)
  
  gset = scdb_gset(atlas$gset_id)
  ref_mc2d = scdb_mc2d(atlas$mc2d_id)
  
  gene_name_map = gen_10x_mars_gene_match(mars_mc_id = mc_id, tenx_mc_id = atlas$mc_id)
  #check cor of all cells, features
  common_genes_ref = names(gset@gene_set)
  common_genes_ref = common_genes_ref[!is.na(gene_name_map[common_genes_ref])]
  common_genes_ref = intersect(common_genes_ref, rownames(ref_mc@e_gc))
  common_genes_ref = common_genes_ref[!is.na(gene_name_map[common_genes_ref])]
  common_genes_ref = common_genes_ref[gene_name_map[common_genes_ref] %in% rownames(query_mc@e_gc)]
  common_genes = gene_name_map[common_genes_ref]
  
  if(mean(!is.null(common_genes)) < 0.5) {
    stop("less than half of the atlas feature genes can be mapped to reference gene names. Probably should provide a name conversion table")
  }
  
  # next calculate average age of reference atlas metacells
  mc_ref_ag = table(ref_mc@mc,mat_ref@cell_metadata[names(ref_mc@mc),"stage"])
  mc_ref_ag = mc_ref_ag[,1:9]
  mc_ref_ag = mc_ref_ag[rowSums(mc_ref_ag)>0,]
  mc_list= as.numeric(rownames(mc_ref_ag))
  
  mc_ref_ag_c = t(t(mc_ref_ag)/colSums(mc_ref_ag))
  mc_ref_ag_cn = mc_ref_ag_c/rowSums(mc_ref_ag_c)
  
  time_points = 6.5 + c(0:8)/4 
  mc_mean_age = (mc_ref_ag_cn %*% time_points)[,1]
  names(mc_mean_age) = rownames(mc_ref_ag_cn)
  
  # next find best reference metacell for each single cell from the query matrix
  feats = mat@mat[common_genes, names(query_mc@mc)]
  rownames(feats) = common_genes_ref
  
  ref_abs_lfp = log(1e-6+ref_mc@e_gc[common_genes_ref,as.character(mc_list)])
  ref_abs_fp = ref_mc@e_gc[common_genes_ref,as.character(mc_list)]
  
  
  cross = tgs_cor((cbind(ref_abs_lfp, as.matrix(feats))))
  cross1 = cross[1:ncol(ref_abs_lfp),ncol(ref_abs_lfp)+1:ncol(feats)]
  
  f = rowSums(is.na(cross1))
  cross1[f,] = 0
  
  best_ref = as.numeric(unlist(apply(cross1,2,function(x) names(which.max(x)))))
  best_ref_cor = apply(cross1, 2, max)
  
  # only include cells in age calculation that are projected on an included cell type from the reference atlas 
  list_of_ignored_reference_mcs = which(ref_mc@colors %in% list_of_ignored_colors)
  
  best_ref_mean_age = mc_mean_age[as.character(best_ref)]
  names(best_ref_mean_age) = colnames(cross1)
  best_ref_mean_age = best_ref_mean_age[!(best_ref %in% list_of_ignored_reference_mcs)]
  
  
  # summarise single embryo time distribution in a data.frame
  emb_age_dist = split(best_ref_mean_age,mat@cell_metadata[names(best_ref_mean_age),"embryo"])
  n_cells = sapply(emb_age_dist,length)
  mean_age = sapply(emb_age_dist,mean)
  median_age = sapply(emb_age_dist,median)
  var_age = sapply(emb_age_dist,var)
  
  sing_emb_age = data.frame(embryo = names(emb_age_dist),
                            ref_gastru_atlas_age = signif(mean_age[names(emb_age_dist)],4),
                            ref_gastru_atlas_rank = rank(mean_age[names(emb_age_dist)]),stringsAsFactors = F)
  
  sing_emb_age  = sing_emb_age[order(sing_emb_age$ref_gastru_atlas_rank),]
  
  
  
  return(sing_emb_age)
}

wt10_gen_intrinsinc_ranking = function() {
  
  mat = scdb_mat("sing_emb_wt10")
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  legc = log2(1e-5 + mc@e_gc)
  
  # define blood mcs
  f = ( legc["Cited4",] > -14 )
  blood_mcs = which(f)
  
  # next extraembryonic endo
  x1 = -8.4
  y1 = -14
  x2 = -11
  y2 = -8.4
  
  b_exe_endo = (y2 - y1)/(x2 - x1)
  a_exe_endo = (y1*x2 - y2*x1)/(x2 - x1)
  
  f = legc["Ttr",] > a_exe_endo + b_exe_endo*legc["Apoe",]
  exe_endo_mcs = which(f)
  
  excluded_mcs = c(blood_mcs,exe_endo_mcs)
  
  df_init_ord = read.table("data/intrinsic_temporal_ranking/intrinsic_ranking_embryo_initial_order_first_iteration.txt",sep = "\t",stringsAsFactors = F)
  
  emb_age_group = unique(mat@cell_metadata[names(mc@mc),c("embryo","age_group")])
  
  
  intrinsic_ranking = sing_emb_intrinsic_ranking(mat_id = "sing_emb_wt10",
                             graph_id = "sing_emb_wt10",
                             mc_id = "sing_emb_wt10_recolored",
                             df_embryo_coarse_ranking = df_init_ord,
                             excluded_mcs = excluded_mcs)
  
  
  intrinsic_ranking = sing_emb_intrinsic_ranking(mat_id = "sing_emb_wt10",
                                                 graph_id = "sing_emb_wt10",
                                                 mc_id = "sing_emb_wt10_recolored",
                                                 df_embryo_coarse_ranking = intrinsic_ranking$embryo_final_order,
                                                 excluded_mcs = excluded_mcs)
  
  return(intrinsic_ranking)
}

sing_emb_intrinsic_ranking = function(mat_id, graph_id,mc_id,df_embryo_coarse_ranking,excluded_mcs,
                                      number_of_neighbors = 50) {
  
  # This function is based only on a coarse initial ranking
  cgraph = scdb_cgraph(graph_id)
  mat = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  number_of_neighbors = number_of_neighbors
  remove_neighbors_from_same_batch = TRUE
  reshuffle = T
  
  # initial ranking - use coarse morphological marks
  
  #check if the embryo names in df_embryo_coarse_ranking are all contained in the mat
  embryos = df_embryo_coarse_ranking$embryo
  embryos = embryos[embryos %in% unique(mat@cell_metadata[colnames(mat@mat),"embryo"])]
  
  included_mcs =  setdiff(c(1:ncol(mc@e_gc)),excluded_mcs)
  included_cells = names(mc@mc)[mc@mc %in% included_mcs]
  
  #embryos = as.character(unique(mat@cell_metadata[names(mc@mc),"embryo"]))
  
  cell_to_embryo = as.character(mat@cell_metadata[names(mc@mc),"embryo"])
  names(cell_to_embryo) = names(mc@mc)
  
  #old way
  #included_edges = cgraph@edges[(( cgraph@edges$mc1 %in% names(mc@mc) ) & ( cgraph@edges$mc2 %in% names(mc@mc) )),]
  #included_edges = included_edges[included_edges$w > (1 - number_of_neighbors/100),]
  #included_edges = included_edges[( included_edges$mc1 %in% included_cells ) & ( included_edges$mc2 %in% included_cells ),]
  
  #new way
  included_edges = cgraph@edges[(( cgraph@edges$mc1 %in% included_cells ) & ( cgraph@edges$mc2 %in% included_cells )),]
  included_edges$mc1 = as.character(included_edges$mc1)
  included_edges$mc2 = as.character(included_edges$mc2)
  
  # exclude embryo self edges
  f = !(cell_to_embryo[included_edges$mc1] == cell_to_embryo[included_edges$mc2])
  included_edges = included_edges[f,]

  included_edges = included_edges[order(included_edges$mc1),]
  list_temp = tapply(X = 1-included_edges$w,INDEX = as.character(included_edges$mc1),FUN = function(x) {rank(x)})
  neighbor_ranks = unlist(tapply(X = 1-included_edges$w,INDEX = as.character(included_edges$mc1),FUN = rank))
  included_edges = included_edges[neighbor_ranks < number_of_neighbors +1,]
  included_cells = unique(included_edges$mc1)
  
  mat_emb_neighbors = table(mat@cell_metadata[as.character(included_edges$mc1),"embryo"],
                            mat@cell_metadata[as.character(included_edges$mc2),"embryo"])
  mat_emb_neighbors = mat_emb_neighbors[embryos,embryos]
  mat_emb_neighbors = as.matrix(mat_emb_neighbors)
  #mat_emb_neighbors = mat_emb_neighbors - diag(diag(mat_emb_neighbors))
  knn_neighbors = rowSums(mat_emb_neighbors)
  batch_to_embryo_table = unique(mat@cell_metadata[names(mc@mc),c("embryo","Sort.Date")])
  
  embryo_weight = table(mat@cell_metadata[included_cells,"embryo"])
  embryo_weight = embryo_weight[embryos]
  
  # Normalize each column by the total number of cells this embryo has.
  mat_emb_neighbors = t(t(mat_emb_neighbors)/as.numeric(embryo_weight))

  
  if (remove_neighbors_from_same_batch) {
    for (i in 1:length(embryos)) {
      embryo = embryos[i]
      batch_id = batch_to_embryo_table$Sort.Date[batch_to_embryo_table$embryo == embryo]
      embryos_same_batch = as.character(batch_to_embryo_table$embryo[batch_to_embryo_table$Sort.Date == batch_id])
      #mat_emb_neighbors[embryos_same_batch,embryo] = 0
      
      embryo_neighbors = setdiff(embryos[max(1,i-10):min(i+10,length(embryos))],embryo)
      # next comes the renormalization part.
      if(length(intersect(embryo_neighbors,embryos_same_batch)) > 0) {
        f = intersect(embryo_neighbors,embryos_same_batch)
        same_batch_mean_of_neighbors = mean(mat_emb_neighbors[embryo,f])
        f2 = setdiff(embryo_neighbors,embryos_same_batch)
        if (length(f2) > 0) {
          other_batches_mean_of_neighbors = mean(mat_emb_neighbors[embryo,f2])
          if (same_batch_mean_of_neighbors > 0) {
            mat_emb_neighbors[embryo,embryos_same_batch] = mat_emb_neighbors[embryo,embryos_same_batch]*
              other_batches_mean_of_neighbors/same_batch_mean_of_neighbors
          }
        }
      }
      
    }
  }
  
  mat_emb_neighbors = mat_emb_neighbors/rowSums(mat_emb_neighbors)
  
  embryo_init_order = embryos
  mat_emb_init = mat_emb_neighbors
  
  # calculate reshuffled matrix
  out = sing_emb_time_ranking_shuffle_rows_and_cols(embryos,mat_emb_neighbors)
  
  mat_emb_neighbors = out$mat_reshuffled
  embryos = out$embryos_final_order
  
  embryo_final_order = data.frame(embryo = rownames(mat_emb_neighbors),transcriptional_rank = c(1:nrow(mat_emb_neighbors)),
                                  stringsAsFactors = F)
  
  
  return(list(similarity_mat = mat_emb_neighbors,embryo_final_order = embryo_final_order))
}


sing_emb_time_ranking_shuffle_rows_and_cols = function(embryos,mat_emb_neighbors) {
  
  count = 0
  n_emb = length(embryos)
  n_iter_max = 100
  n_iter = 0
  delta_count = 1
  
  while ((n_iter < n_iter_max)  & (delta_count > 0) ) {
    
    old_count = count
    n_iter = n_iter + 1
    
    for(i in 1:(n_emb-1)) {
      
      embryo = embryos[i]
      next_embryo = embryos[i+1]
      
      if (i == 1) {
        row_2 = sum(mat_emb_neighbors[embryo,(i+2):n_emb]) - sum(mat_emb_neighbors[next_embryo,(i+2):n_emb])
        col_2 = sum(mat_emb_neighbors[(i+2):n_emb,embryo]) - sum(mat_emb_neighbors[(i+2):n_emb,next_embryo])
        if (row_2 +col_2 > 0) {
          embryos[i] = next_embryo
          embryos[i+1] = embryo
          mat_emb_neighbors = mat_emb_neighbors[embryos,embryos]
          count = count + 1 
        }
      } else if (i == (n_emb-1)) {
        row_1 = sum(mat_emb_neighbors[embryo,1:(i-1)]) - sum(mat_emb_neighbors[next_embryo,1:(i-1)])
        col_1 = sum(mat_emb_neighbors[1:(i-1),embryo]) - sum(mat_emb_neighbors[1:(i-1),next_embryo])
        if (row_1 + col_1 > 0) {
          embryos[i] = next_embryo
          embryos[i+1] = embryo
          mat_emb_neighbors = mat_emb_neighbors[embryos,embryos]
          count = count + 1 
        }
        
        
      } else {
        row_1 = sum(mat_emb_neighbors[embryo,1:(i-1)]) - sum(mat_emb_neighbors[next_embryo,1:(i-1)])
        row_2 = sum(mat_emb_neighbors[embryo,(i+2):n_emb]) - sum(mat_emb_neighbors[next_embryo,(i+2):n_emb])
        col_1 = sum(mat_emb_neighbors[1:(i-1),embryo]) - sum(mat_emb_neighbors[1:(i-1),next_embryo])
        col_2 = sum(mat_emb_neighbors[(i+2):n_emb,embryo]) - sum(mat_emb_neighbors[(i+2):n_emb,next_embryo])
        if (row_2 - row_1 + col_2 - col_1 > 0) { 
          embryos[i] = next_embryo
          embryos[i+1] = embryo
          mat_emb_neighbors = mat_emb_neighbors[embryos,embryos]
          count = count + 1 
        }
      }
      
      
    }
    delta_count = count - old_count
    print(delta_count)
  }
  
  
  return(list(mat_reshuffled = mat_emb_neighbors,embryos_final_order = embryos))
} 

