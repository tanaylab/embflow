# script for downloading and initializing the single-cell RNA-seq database for the metacell package

download.file("https://embflow.s3.eu-west-1.amazonaws.com/data_embflow.tar.gz", "data_embflow.tar.gz")

system("tar -xvf data_embflow.tar.gz data")

download.file("https://embflow.s3.eu-west-1.amazonaws.com/scrna_db_embflow.tar.gz","scrna_db_embflow.tar.gz")

system("tar -xvf scrna_db_embflow.tar.gz scrna_db")

file.remove("data_embflow.tar.gz")
file.remove("scrna_db_embflow.tar.gz")


if(!dir.exists("figs")) {
  dir.create("figs")
} 
if(!dir.exists("figs/paper_figs")) {
  dir.create("figs/paper_figs")
}

