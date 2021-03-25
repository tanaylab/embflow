# script for downloading and initializing the single-cell RNA-seq database for the metacell package

download.file("http://www.wisdom.weizmann.ac.il/~atanay/embflow/data_embflow.tar.gz", "data_embflow.tar.gz")

system("tar -xvf data_embflow.tar.gz data")

download.file("http://www.wisdom.weizmann.ac.il/~atanay/embflow/scrna_db_embflow.tar.gz","scrna_db_embflow.tar.gz")

system("tar -xvf scrna_db_embflow.tar.gz scrna_db")

file.remove("data_embflow.tar.gz")
file.remove("scrna_db_embflow.tar.gz")


if(!dir.exists("figs")) {
  dir.create("figs")
} 
if(!dir.exists("figs/paper_figs")) {
  dir.create("figs/paper_figs")
}

