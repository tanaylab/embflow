# summary of preprocessing of Foxc Control plates


# 1. cells with #UMIs < 1000 were removed
# 2. the following genes were removed from the single cell matrix
# mat = scdb_mat(mat_name)
# nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
# bad_genes = c(grep("^mt\\-", nms, v=T), "Neat1",grep("ERCC", nms,v=T), "Atpase6", "Xist", "Malat1", "Cytb","AK018753","AK140265","AK163440","DQ539915")
# mcell_mat_ignore_genes(mat_name, mat_name,bad_genes,reverse=F) 
# mcell_mat_ignore_small_cells(mat_name, mat_name, 1000)

# 3. duplicate cells (from FACS sorting) were removed
# 4. cells from empty wells were removed (should anyway have less than 1000 UMIs)

# 5. cells were gated based on FACS GFP fluorescent signal using the script foxc_control_gating.r
# source("scripts/foxc12/preprocessing/foxc_control_gating.r")
# gating_of_foxc_control("figs/paper_figs/fig_s6/fig_s6b")

# 6. Extraembryonic ectoderm and parietal endoderm cells were removed from analysis
# see code in scripts/foxc12/preprocessing/foxc_control_remove_exe_ectoderm_and_parietal_endo_cls.r

# 7. Single cell matrix was split into four submatrices, one for Foxc/Control Chimera/Tetraploid
#    Each matrix was merged with the wildtype sing_emb_wt10 matrix
#    see  scripts/foxc12/preprocessing/foxc_control_split_into_experiment_type.r