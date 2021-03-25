library("metacell")
scdb_init("scrna_db/",force_reinit = T)

# permutation test for 2 tfs


# for every feat gene compute the two most predictive TFs 
# among them, select the one that has the highest single predictive power and fix it
# next randomize the rows of tfs * metacell matrix excluding the selected TFS
# find best predictive pair from randomized data. 
# repeat a 500 times
# create boxplot of randomized R squared vs real R squared

gen_fig_s5_plots = function() {
  
  if(!dir.exists("figs/paper_figs/fig_s5")) {
    dir.create("figs/paper_figs/fig_s5")
  }
  
  fig_s5a()
  # If the following files 
  #
  # r_sq_permut_test_tmp_1.Rda
  # r_sq_permut_test_tmp_2.Rda
  # r_sq_permut_test_tmp_3.Rda
  #
  # are not available in the subfolder data/fig5 please rerun
  # permutation_test_two_tf_pred() with the input values 1,2,3:
  #
  # permutation_test_two_tf_pred(1)
  # permutation_test_two_tf_pred(2)
  # permutation_test_two_tf_pred(3)
  
  fig_s5b()
  fig_s5c()
}

fig_s5b = function() {
  
  best_predict = read.table(file = "data/fig5/best_pred_2tfs.txt",sep = "\t",stringsAsFactors = F)
  best_predict_1tf = read.table(file = "data/fig5/best_pred_1_tf.txt",sep = "\t",stringsAsFactors = F)
  colnames(best_predict_1tf) = c("targ","best_tf_1tf","rsq_one_tf")
  best_predict = left_join(best_predict,best_predict_1tf, by = "targ")
  
  load("data/fig5/r_sq_permut_test_tmp_1.Rda")
  r_sq_ls_all = list()
  r_sq_ls_all[c(1:70)] = r_sq_ls[1:70]
  load("data/fig5/r_sq_permut_test_tmp_2.Rda")
  r_sq_ls_all[c(71:140)] = r_sq_ls[71:140]
  load("data/fig5/r_sq_permut_test_tmp_3.Rda")
  r_sq_ls_all[c(141:nrow(best_predict))] = r_sq_ls[141:nrow(best_predict)]
  
  names(r_sq_ls_all) = best_predict$targ
  
  delta_r_sq = best_predict$m_rsq - best_predict$rsq_one_tf
  names(delta_r_sq) = best_predict$targ
  delta_r_sq = delta_r_sq[order(delta_r_sq,decreasing = T)]
  r_sq_ls_all = r_sq_ls_all[names(delta_r_sq)]
  

  names(r_sq_ls_all) = substr(names(r_sq_ls_all),1,10)
  
  
  ylim_box = c(0,max(unlist(r_sq_ls_all),delta_r_sq))
  n_min = 1
  n_max = 100
  pdf("figs/paper_figs/fig_s5/fig_s5b_top100.pdf",useDingbats = F,w = 20, h = 5)
  par(mar = c(6,6,4,4))
  boxplot(r_sq_ls_all[n_min:n_max],ylim = ylim_box,las = 2,ylab = expression(paste(Delta,paste(" R"^"2"))),main = expression(paste("Difference in ",paste("R"^"2")," values using 2 TFs vs 1 TF in linear model regression")))
  points(delta_r_sq[n_min:n_max],col = "red",pch = 19)
  dev.off()
  
  n_min = 101
  n_max = nrow(best_predict)
  pdf("figs/paper_figs/fig_s5/fig_s5b_top200.pdf",useDingbats = F,w = 20, h = 5)
  par(mar = c(6,6,4,4))
  boxplot(r_sq_ls_all[n_min:n_max],ylim = ylim_box,las = 2,ylab = expression(paste(Delta,paste(" R"^"2"))))
  points(delta_r_sq[n_min:n_max],col = "red",pch = 19)
  dev.off()
  
  
  
  
}


permutation_test_two_tf_pred = function(n_ind) {
  

  
  best_fit_1_tf_1_fixed = function(targ,tf_list,tf_fixed,legc,legc_shuff,best_one_tf) {
    r_sq = 0
    best_fit = c(0,0,0)
    tf_list_f = setdiff(tf_list,tf_fixed)
    for(i in 1:length(tf_list_f)) {
      tf1 = tf_list[i]
      m = summary(lm(legc[targ,] ~ legc_shuff[tf1,] + legc[tf_fixed,]))
      if(m$r.squared > r_sq) {
        best_fit = c(targ,tf1, m$r.squared)
        r_sq = m$r.squared
      }
      m2 = summary(lm(legc[targ,] ~ legc_shuff[tf1,] + legc[best_one_tf,]))
      if(m2$r.squared > r_sq) {
        best_fit = c(targ,tf1, m2$r.squared)
        r_sq = m2$r.squared
      }
    }
    return(r_sq)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(mc@e_gc[,mc_f] + 3e-5)
  
  best_predict = read.table(file = "data/fig5/best_pred_2tfs.txt",sep = "\t",stringsAsFactors = F)
  best_predict_1tf = read.table(file = "data/fig5/best_pred_1_tf.txt",sep = "\t",stringsAsFactors = F)
  colnames(best_predict_1tf) = c("targ","best_tf_1tf","rsq_one_tf")
  best_predict = left_join(best_predict,best_predict_1tf, by = "targ")
  
  best_predict$rsq_tf_fix =  ifelse(best_predict$tf1_r_sq > best_predict$tf2_r_sq,best_predict$tf1_r_sq,best_predict$tf2_r_sq)
  
  markers = best_predict$targ
  
  tf_list = read.table("data/fig5/tf_list_mesoderm.txt",sep = "\t",stringsAsFactors = F)$x
  
  set.seed(112)
  r_sq_ls = list()
  
  if (n_ind == 1) {
    ind_ls = c(1:70)
  } else if (n_ind == 2) {
    ind_ls = c(71:140)
  } else {
    ind_ls = c(141:nrow(best_predict))
  }
  
  for (i in ind_ls) {
    print(i)
    targ = best_predict$targ[i]
    r_sq_tf1 = best_predict$tf1_r_sq[i]
    r_sq_tf2 = best_predict$tf2_r_sq[i]
    
    tf_fix = best_predict$tf1[i]
    if(r_sq_tf2 > r_sq_tf1) {
      tf_fix = best_predict$tf2[i]
    }
    rsq_tf_fix = best_predict$rsq_tf_fix[i]
    best_one_tf = best_predict$best_tf_1tf[i]
    
    
    r_sq_v = sapply(c(1:100),function(n) {
      
      mc_shuff = sample(ncol(legc))
      legc_shuff = legc[tf_list,mc_shuff]
      
      r_sq = best_fit_1_tf_1_fixed(targ = targ,tf_list = tf_list,tf_fixed = tf_fix,legc = legc,legc_shuff = legc_shuff,best_one_tf = best_one_tf)
      return(r_sq)
    })
    
    r_sq_ls[[i]] = r_sq_v - best_predict$rsq_one_tf[i]
    
  }
  
  delta_r_sq = best_predict$m_rsq[ind_ls] - best_predict$rsq_one_tf[ind_ls]
  
  
  save(r_sq_ls,file = sprintf("data/fig5/r_sq_permut_test_tmp_%d.Rda",n_ind))
  
}




fig_s5a = function() {
  
  tf_list = read.table("data/fig5/tf_list_mesoderm.txt",sep = "\t",stringsAsFactors = F)$x
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(1e-5+mc@e_gc[tf_list,mc_f])
  
  tf_max = apply(legc,1,max) - log2(1e-5)
  
  tf_max = sort(tf_max,decreasing = T)
  
  pdf("figs/paper_figs/fig_s5/fig_s5a.pdf",useDingbats = F,width = 12,h = 4.5)
  barplot(tf_max,las = 2,yaxt = 'n',ylim =  c(0,-8-log2(1e-5)))
  axis(side = 2,labels = c(-16.6,-14,-12,-10,-8),at = c(-16.6,-14,-12,-10,-8) - log2(1e-5))
  dev.off()
  
}

fig_s5c = function(plot_pdf = T) {
  
  fig_dir = "figs/paper_figs/fig_s5/fig_s5c"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac mesoderm","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(mc@e_gc[,mc_f] + 3e-5)
  
  # these external tables are generated by the function predict_variable_genes_lm()
  best_pred_1_tf = read.table(file = "data/fig5/best_pred_1_tf.txt",sep = "\t",stringsAsFactors = F)
  best_pred_2_tfs = read.table(file = "data/fig5/best_pred_2tfs.txt",sep = "\t",stringsAsFactors = F)
  
  rownames(best_pred_1_tf) = best_pred_1_tf$targ
  rownames(best_pred_2_tfs) = best_pred_2_tfs$targ
  
  target_genes = c("Col23a1","Dll3","Hsd11b2","Tmem88","Vim")
  
  for (targ in target_genes) {
    
    tf1 = best_pred_2_tfs[targ,"tf1"]
    tf2= best_pred_2_tfs[targ,"tf2"]
    
    if(plot_pdf) {
      pdf(sprintf("%s/%s_vs_%s.pdf",fig_dir,targ,tf1))
    } else {
      png(sprintf("%s/%s_vs_%s.png",fig_dir,targ,tf1),w = 500,h = 500)
    }
    par(mar = c(6,6,2,4))
    plot(x = legc[tf1,],y = legc[targ,],pch = 19,cex = 4,col = mc@colors[mc_f],xlab = tf1,ylab = targ,cex.lab = 2)
    dev.off()
    
    if (plot_pdf) {
      pdf(sprintf("%s/%s_vs_%s.pdf",fig_dir,targ,tf2))
    } else {
      png(sprintf("%s/%s_vs_%s.png",fig_dir,targ,tf2),w = 500,h = 500)
    }
    par(mar = c(6,6,2,4))
    plot(x = legc[tf2,],y = legc[targ,],pch = 19,cex = 4,col = mc@colors[mc_f],xlab = tf2,ylab = targ,cex.lab = 2)
    dev.off()
    
    out = lm(legc[targ,] ~ legc[tf1,]+legc[tf2,])
    
    if(plot_pdf) {
      pdf(sprintf("%s/%s_vs_%s_and_%s.pdf",fig_dir,targ,tf1,tf2))
    } else {
      png(sprintf("%s/%s_vs_%s_and_%s.png",fig_dir,targ,tf1,tf2),w = 500,h = 500)
    }
    par(mar = c(6,6,2,4))
    plot(x = out$fitted.values,y = legc[targ,],pch = 19,cex = 4,col = mc@colors[mc_f],xlab = sprintf("Prediction from %s and %s",tf1,tf2),ylab = targ,cex.lab = 2)
    dev.off()
    
  }
  
  
}

