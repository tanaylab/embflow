A single embryo, single cell time-resolved model for mouse gastrulation
=======================================================================

This repository contains all the code for reproducing the analysis from the gastrulation flow paper Mittnenzweig et al. (2021). The core analysis is done with the [metacell](https://github.com/tanaylab/metacell) R package, that also contains the code for generating the network flow model.

### Quick links

- Metacell paper: Baran et al. 2019 [Genome Biol](https://doi.org/10.1186/s13059-019-1812-2)

- [Metacell](https://github.com/tanaylab/metacell) R package

- Raw FASTQ files and processed UMI tables are available under GEO accession [GSE169210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169210)

### Requirements

- metacell package
- lpsymphony
- tidyverse
- pheatmap
- gridExtra
- Matrix
- tidyverse
- shape
- umap

### Usage

After cloning the github repository, users should open an R session in the repository root directory and download/initialize the scRNA database (~ 4.7 GB): 
``` r
# Loading code and downloading required data files
source("scripts/init_db.r")

```
The repository root directory should now contain the subfolders *scripts/* containing all the R scripts, *scrna_db* containing the metacell R objects, *config/*, *data/* containing additional data generated by the scripts and *figs/paper_figs/*.

#### Regenerating plots for a specific figure
For each figure (Figures 1-7 and S1-7), there is a corresponding script in *scripts/generate_paper_figures/*. Each script contains a function *gen_fig_xyz_plots()* at the top, that contains further subfunctions and explanations related to the analysis of that figure. E.g., for regenerating the plots of figure 1, users should run the following code:

``` r
# load metacell package
library("metacell")
# initializing the metacell scrna database
scdb_init("scrna_db")

# Generating plots of Figure 1
source("scripts/generate_paper_figures/fig_1.r")
gen_fig_1_plots()

```
The content of *gen_fig_1_plots()* looks as follwos:
``` r
gen_fig_1_plots = function() {
  
  if(!dir.exists("figs/paper_figs")) {
    dir.create("figs/paper_figs")
  }
  dir_name = "figs/paper_figs/fig1"
  
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  fig1_b()
  fig1_cde()
  fig1_f()
  fig1_g_mc_time_distributions()
  fig1_g_heatmap(plot_pdf = T)
  fig1_h()
  
}
```
Figure plots are saved in *figs/paper_figs/fig1/*


#### Computing metacell object, manifold graph and network flow model for wildtype manifold
Standard metacell analysis is performed as described in [Baran et al. 2019](https://doi.org/10.1186/s13059-019-1812-2). To recompute the metacell object, please run
``` r
source("scripts/generate_mc_mgraph_network/gen_mc.r")
```
This will generate a metacell object with id *sing_emb_wt10_bs500f*. Note, that because of random seeding of the boostrap procedure involved in calculating the metacell cover,
the computed metacell cover will slightly deviate from *sing_emb_wt10_recolored* used in the paper. Manifold graphs and 2D projections can be recomputed through
``` r
source("scripts/generate_mc_mgraph_network/gen_mgraph.r")
source("scripts/generate_mc_mgraph_network/gen_mgraph_umap.r")
```
The network flow model can be generated using
``` r
source("scripts/generate_mc_mgraph_network/gen_network.r")
build_sing_emb_wt10_network()
```
Metacells were clustered and annotated using the network flow model.
``` r
source("scripts/generate_mc_mgraph_network/annot_mc_by_flows.r")
cluster_metacells_by_flow(mct_id = "sing_emb_wt10",K = 65)
```

#### Single-embryo timing
To regenerate the single-embryo timing data underlying Figure 1, please run
``` r
# load metacell package
library("metacell")
# initializing the metacell scrna database
scdb_init("scrna_db")

source("scripts/single_embryo_timing.r")
gen_fig_1_plots()
embryo_ranks = gen_single_embryo_timing()

# subfunctions calculating intrinsic_rank and reference_rank of each embryo
# are contained in gen_single_embryo_timing()
```
The output data frame *embryo_ranks* was added to the single-cell metadata information of the metacell matrix object. All subsequent functions using single-embryo time information, are extracting it from the *cell_metadata* entry of the WT metacell single-cell matrix object *sing_emb_wt10*.
``` r
# load metacell package
library("metacell")
# initializing the metacell scrna database
scdb_init("scrna_db")

mat = scdb_mat("sing_emb_wt10")
md = mat@cell_metadata
```

#### Parameter stability analysis of network flow model
The parameter stability analysis of network flows underlying Figure S2A can be regenerated using
``` r
# regnerate data - this might take some time
source("scripts/parameter_stability_analysis.r")
gen_parameter_stability_analysis()

# replotting Figure S2A
source("scripts/generate_paper_figures/fig_s2.r")
fig_s2a()
```

#### Foxc12 chimera and tetraploid analysis
To generate specific plots of Figures 6, S6 and S7, please run the corresponding functions from *fig_6.r*, *fig_s6.r* or *fig_s7.r*. Users interested in recomputing parts of the Foxc12 chimera and tetraploid embryo analysis (not needed for regenerating the plots), should run the following functions:
``` r
library("metacell")
scdb_init("scrna_db/")

source("scripts/foxc12/generate_chimera_tetraploid_data_analysis.r")

# Chimera embryos injected with Foxc12 DKO cells
foxc_chimera_generate_time_and_cell_type_annotation()

# Chimera embryos injected with control cells
control_chimera_generate_time_and_cell_type_annotation()

# Tetraploid embryos injected with Foxc12 DKO cells
foxc_tetraploid_generate_time_and_cell_type_annotation()

# Tetraploid embryos injected with control cells
control_tetraploid_generate_time_and_cell_type_annotation()
```
This will transfer cell-type and time annotation from the wt atlas to chimera/tetraploid embryos. Output is saved in *data/chimera_tetraploid_analysis/*. Scripts involved in preprocessing plates from the chimera and tetraploid embryo analyis are saved in the *scripts/foxc12/preprocessing/*. This includes
- Gating of single cells using the FACS GFP channel
- Removing cells from extraembryonic ectoderm and parietal endoderm
- Merging each single-cell matrix with the wt single-cell matrix and creating a joint single-cell graph (metacell *cgraph* object).
See summary_preprocessing.r and the corresponding scripts for more details.


### Metacell R objects
- 
-
- 
- 
- 
-
- 
- 
- 
- 
- 
- 
- 
-
- 
-
-
- 


### Contact
For help, please contact <markus.mittnenzweig@weizmann.ac.il>