library("umap")
tgconfig::override_params(config_file = "config/sing_emb_mgraph.yaml",package = "metacell")
# generate mgraph and mc2d object

# generate mgraph using logistic distance
# In config/sing_emb_mgraph.yaml use parameters
# mcell_mgraph_T_edge: 0.001
# mcell_mgraph_max_confu_deg: 4

feat_gset = "sing_emb_wt10"
mc_id = "sing_emb_wt10_recolored"
mgraph_id = paste0(mc_id,"_logist")
logist_loc = 1
logist_scale = 0.2
logist_eps = 4e-5
max_d_fold = 3

mc = scdb_mc(mc_id)
gset = scdb_gset(feat_gset)
feat_genes = names(gset@gene_set)

mgraph = mc2d_comp_mgraph_param(mc, feat_genes, logist_loc, logist_scale, logist_eps, max_d_fold)

scdb_add_mgraph(id = mgraph_id,mgraph = tgMCManifGraph(mc_id = mc_id,mgraph = mgraph))


# next generate 2d projection using umap

tgconfig::set_param(param = "mcell_mc2d_max_confu_deg",value = 4,package = "metacell")
mgraph = mc2d_comp_mgraph_param(mc, feat_genes, logist_loc, logist_scale, logist_eps, max_d_fold)

graph_id = "sing_emb_wt10"
mc2d_id = paste0(mc_id,"_umap")
symmetrize = F
umap_mgraph = F

uconf = umap.defaults
#uconf$n_neighbors=6
#uconf$min_dist=0.9
uconf$n_neighbors=4
uconf$min_dist =0.9
uconf$bandwidth=1.3


#mgraph = scdb_mgraph(mgraph_id)
#mgraph = mgraph@mgraph

mc_xy = mc2d_comp_graph_coord_umap(mc, feat_genes, mgraph, uconf, umap_mgraph)
xy = mc2d_comp_cell_coord(mc_id = mc_id,graph_id =  graph_id, mgraph = mgraph, cl_xy = mc_xy, symmetrize=symmetrize)
scdb_add_mc2d(mc2d_id, tgMC2D(mc_id, mc_xy$mc_x, mc_xy$mc_y, xy$x, xy$y, mgraph))

mcell_mc2d_plot(mc2d_id)
