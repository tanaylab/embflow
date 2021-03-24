# definition of cell types that belong to endoderm, ectoderm or embryonic mesoderm respectively

embryonic_meso_ct_colors = function() {
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  ct_to_col = mc_wt@color_key$color
  names(ct_to_col) = mc_wt@color_key$group
  
  included_ct = c("Caudal mesoderm","Early nascent mesoderm","Late nascent mesoderm","Paraxial mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm",
                  "Cardiac mesoderm","Cardiomyocytes")
  
  if(sum(is.na(ct_to_col[included_ct])) > 0) {
    stop(sprintf("one of the cell types is not correctly defined \n %s",included_ct[is.na(ct_to_col[included_ct])]))
  }
  
  return(ct_to_col[included_ct])
}

endo_ct_colors = function() {
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  ct_to_col = mc_wt@color_key$color
  names(ct_to_col) = mc_wt@color_key$group
  
  included_ct = c("Definitive endoderm","Foregut","Node/Notochord")
  
  
  if(sum(is.na(ct_to_col[included_ct])) > 0) {
    stop("one of the cell types is not correctly defined")
  }
  
  return(ct_to_col[included_ct])
}

ectoderm_ct_colors = function() {
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  ct_to_col = mc_wt@color_key$color
  names(ct_to_col) = mc_wt@color_key$group
  
  included_ct = c("Surface ectoderm","Rostral neural plate","Definitive ectoderm","Caudal neural plate","Caudal epiblast")
  
  if(sum(is.na(ct_to_col[included_ct])) > 0) {
    stop("one of the cell types is not correctly defined")
  }
  
  return(ct_to_col[included_ct])
}


