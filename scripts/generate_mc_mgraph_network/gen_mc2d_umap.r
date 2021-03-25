library("umap")
# generate mgraph and mc2d object

gen_mc2d_umap_wt10 = function() {
  
  mc_id = "sing_emb_wt10_recolored"
  mgraph = scdb_mgraph("sing_emb_wt10_recolored_logist")
  mc = scdb_mc(mc_id)
  gset = scdb_gset("sing_emb_wt10")
  graph_id = "sing_emb_wt10"
  feat_genes = names(gset@gene_set)
  # next generate 2d projection using umap
  
  tgconfig::set_param(param = "mcell_mc2d_max_confu_deg",value = 4,package = "metacell")
  
  mc2d_id = "sing_emb_wt10_recolored_umap"
  symmetrize = F
  umap_mgraph = F
  
  uconf = umap.defaults
  #uconf$n_neighbors=6
  #uconf$min_dist=0.9
  uconf$n_neighbors=4
  uconf$min_dist =0.9
  uconf$bandwidth=1.3
  
  mc_xy = mc2d_comp_graph_coord_umap(mc, feat_genes, mgraph@mgraph, uconf, umap_mgraph)
  xy = mc2d_comp_cell_coord(mc_id = mc_id,graph_id =  graph_id, mgraph = mgraph@mgraph, cl_xy = mc_xy, symmetrize=symmetrize)
  scdb_add_mc2d(mc2d_id, tgMC2D(mc_id, mc_xy$mc_x, mc_xy$mc_y, xy$x, xy$y, mgraph@mgraph))
  
  mcell_mc2d_plot(mc2d_id)
}




