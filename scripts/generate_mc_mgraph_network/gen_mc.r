library("devtools")
load_all("metacell")
source("paper_scripts/generate_mc_mgraph_network/generic_mc.r")
tgconfig::override_params(config_file = "config/sing_emb.yaml",package = "metacell")

scdb_init("scrna_db",force_reinit = T)
scfigs_init("figs")

# first iteration without out filtered genes
# generate_mc(mat_nm, color_key=NA,recompute = T)
# then run embdyn_find_bad_genes.r and filter bad gene modules
# included bad genes in second iteration
bad_genes = read.table("data/sing_emb_wt10.bad_genes.txt",sep = "\t",stringsAsFactors = F)
bad_genes = bad_genes[,1]

mat_nm = "sing_emb_wt10"
mc_id = paste(mat_nm,"_bs500f",sep="")

generate_mc(mat_nm, color_key=NA,add_bad_genes = bad_genes,recompute = T)
