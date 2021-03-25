library("metacell")
source("scripts/generate_mc_mgraph_network/generic_mc.r")
tgconfig::override_params(config_file = "config/sing_emb.yaml",package = "metacell")

scdb_init("scrna_db",force_reinit = T)
scfigs_init("figs")

# first iteration without out filtered genes
# generate_mc(mat_nm, color_key=NA,recompute = T)
# then remove genes from list of feature genes by filter bad gene modules
# remove bad genes in second iteration
bad_genes = read.table("data/external_data/sing_emb_wt10.bad_genes.txt",sep = "\t",stringsAsFactors = F)
bad_genes = bad_genes[,1]

mat_nm = "sing_emb_wt10"

generate_mc(mat_nm, color_key=NA,add_bad_genes = bad_genes,recompute = T)
