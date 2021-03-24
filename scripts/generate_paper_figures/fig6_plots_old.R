# Fig. 6 plots for the Foxc12 chimera

library("devtools")
load_all("metacell")
tgconfig::override_params(config_file = "config/sing_emb.yaml",package = "metacell")
scdb_init("scrna_db/",force_reinit = T)
scfigs_init("figs")

source("paper_scripts/foxc12/chimera_time_match.R")
source("paper_scripts/foxc12/transfer_color_chimera.r")
source("paper_scripts/foxc12/chim_compare_frequencies.R")
source("paper_scripts/foxc12/diff_expression_analysis.r")


gen_fig6_plots() {
  # first tranfer color and save it
  transfer_color()
  
  # next generate chimera_time_match
  foxc12_wt10_f_time_match_chimeras()
  
  plot_time_distributions(plot_pdf = T)
  
  # plot differences in cell type frequencies
  wt10_foxc12_f_cell_type_frequency_comparison(plot_pdf = T)
  
  # differential expression analysis
  chim_diff_expr()
  diff_expression_per_emb_and_ct_heatmap(included_clusters = c(6,8,9),plot_pdf = T)
  
  chim_plot_single_genes(plot_pdf = T)
}


