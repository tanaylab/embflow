source("scripts/foxc12/transfer_time_annotation.r")
library("Matrix")
library("dplyr")



control_chim_timing = function(tag = "all",included_colors = NULL) {
  
  mat_nm = "control_chim_wt10"
  mc_id = "control_chim_wt10_recolored"
  cgraph_id =  mat_nm
  mat = scdb_mat(mat_nm)
  mc = scdb_mc(mc_id)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")

  excluded_colors = c("#F6BFCB","#7F6874")
  
  if(is.null(included_colors)) {
    included_colors = setdiff(mc_wt@color_key$color,excluded_colors)
  }
  
  load(sprintf("data/chimera_tetraploid_analysis/%s/color_annotation/cmp_annot.Rda",mat_nm))
  query_cls_col = cmp_annot$query_cls_col
  
  f = mat@cell_metadata[colnames(mat@mat),"cell_type"] %in% c("host","control","KO")
  chim_cls = colnames(mat@mat)[f]
  chim_cls = intersect(chim_cls,names(cmp_annot$query_cls_col))
  
  excluded_cls = chim_cls[cmp_annot$query_cls_col[chim_cls] %in% excluded_colors]
  chim_cls = setdiff(chim_cls,excluded_cls)
  
  chim_cls = chim_cls[query_cls_col[chim_cls] %in% included_colors]
  
  chim_embryos = unique(mat@cell_metadata[chim_cls,"embryo"])
   
  wt10_cls = intersect(names(mc_wt@mc)[ !(mc_wt@colors[mc_wt@mc] %in% excluded_colors) ],colnames(mat@mat))
  atlas_time = mat@cell_metadata[wt10_cls,"transcriptional_rank"]
  names(atlas_time) = wt10_cls
 
  data_dir = sprintf("data/chimera_tetraploid_analysis/%s",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  
  data_dir = sprintf("data/chimera_tetraploid_analysis/%s/time_match",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }

  # first get atlas time distribution
  
  atlas_time_dist = get_atlas_time_dist(atlas_time = atlas_time,graph_id = cgraph_id)
  
  # first timing using host cells only
  host_cls = chim_cls[(mat@cell_metadata[chim_cls,"cell_type"] == "host") & ( mat@cell_metadata[chim_cls,"embryo"] %in% chim_embryos )]
  query_cls_md = mat@cell_metadata[host_cls,"embryo"]
  names(query_cls_md) = host_cls
  
  query_time_dist_host = get_query_time_dist(query_cls_md = query_cls_md,atlas_time = atlas_time,graph_id = cgraph_id)
  
  # timing using control cells
  
  control_cls = chim_cls[( mat@cell_metadata[chim_cls,"cell_type"] == "control" ) & ( mat@cell_metadata[chim_cls,"embryo"] %in% chim_embryos )]
  query_cls_md = mat@cell_metadata[control_cls,"embryo"]
  names(query_cls_md) = control_cls
  
  query_time_dist_control = get_query_time_dist(query_cls_md = query_cls_md,atlas_time = atlas_time,graph_id = cgraph_id)
  
  
  time_dist_host_control = list(atlas_time_dist = atlas_time_dist$atlas_time_dist,
                           host = query_time_dist_host$query_time_dist,
                           control = query_time_dist_control$query_time_dist)
  
  if(tag != "all") {
    file_nm =  sprintf("%s/time_dist_host_control_%s.Rda",data_dir,tag)
  } else {
    file_nm =  sprintf("%s/time_dist_host_control.Rda",data_dir)
  }
  
  save(time_dist_host_control,file = file_nm)
  
  chim_emb_summary = as.data.frame.matrix(table(mat@cell_metadata[chim_cls,"embryo"],mat@cell_metadata[chim_cls,"cell_type"]))
  chim_emb_summary$embryo = rownames(chim_emb_summary)
  chim_emb_summary$best_rank_host = NA
  chim_emb_summary[rownames(time_dist_host_control$host),"best_rank_host"] = time_dist_best_match(atlas_time_dist = time_dist_host_control$atlas_time_dist,
                                                                                             query_time_dist = time_dist_host_control$host)
  chim_emb_summary$best_rank_control = NA
  chim_emb_summary[rownames(time_dist_host_control$control),"best_rank_control"] = time_dist_best_match(atlas_time_dist = time_dist_host_control$atlas_time_dist,
                                                         query_time_dist = time_dist_host_control$control)
  
  if(tag != "all") {
    file_nm =  sprintf("%s/time_match_summary_%s.txt",data_dir,tag)
  } else {
    file_nm =  sprintf("%s/time_match_summary.txt",data_dir)
  }
  
  write.table(chim_emb_summary,file = file_nm,sep ="\t",row.names = F)
}

time_dist_best_match = function(query_time_dist,atlas_time_dist) {
  
  query_ref_cor = tgs_cor(t(as.matrix(as.data.frame.matrix(query_time_dist))),
                          t(as.matrix(as.data.frame.matrix(atlas_time_dist))))
  
  best_fit = apply(query_ref_cor,1,FUN = function(x) {
    return(mean(as.numeric(colnames(query_ref_cor))[order(x,decreasing = T)][1:5]))
  })
  
  return(best_fit)
}

control_chimera_plot_cumulative_distribution_control_host = function(mat_nm,fig_dir = NULL,tag = "all",plot_pdf = T,chim_col = "cornflowerblue") {
  
  fig_scale = 1.4
  lwd = 15
  x_pos_text = 85
  xlim_min = 75
  xlim_max = 153
  
  data_dir = sprintf("data/chimera_tetraploid_analysis/%s/time_match",mat_nm)
  if(tag == "all") {
    load(sprintf("%s/time_dist_host_control.Rda",data_dir))
  } else {
    load(sprintf("%s/time_dist_host_control_%s.Rda",data_dir,tag))
  }
  
  host_dist_all = time_dist_host_control$host
  control_dist_all = time_dist_host_control$control
  
  if(tag == "all") {
    fn = sprintf("%s/time_match_summary.txt",data_dir)
  } else {
    fn = sprintf("%s/time_match_summary_%s.txt",data_dir,tag)
  }
  chim_time_summary = read.table(file = fn,sep ="\t",h = T,stringsAsFactors = F)
  chim_embryos = chim_time_summary$embryo[order(chim_time_summary$best_rank_host)]
  
  ks_statistic = rep(0,length(chim_embryos))
  p_value = rep(0,length(chim_embryos))
  delta_median = p_value
  
  for (i in 1:length(chim_embryos)) {
    chimera = chim_embryos[i]
    host_dens = cumsum(host_dist_all[chimera,])/sum(host_dist_all[chimera,])
    control_dens = cumsum(control_dist_all[chimera,])/sum(control_dist_all[chimera,])
    
    host_median = median(rep(c(1:153),host_dist_all[chimera,]))
    control_median = median(rep(c(1:153),control_dist_all[chimera,]))
    
    delta_median[i] = host_median - control_median
    
    ks_control_host_chim = ks.test(x = rep(c(1:153),host_dist_all[chimera,]),y = rep(c(1:153),control_dist_all[chimera,]))
    ks_statistic[i] = ks_control_host_chim$statistic
    p_value[i] = ks_control_host_chim$p.value
    
    
    N_host = sum(host_dist_all[chimera,])
    N_control = sum(control_dist_all[chimera,])
    
    if(!is.null(fig_dir)) {
      if(!dir.exists(fig_dir)) {
        dir.create(fig_dir)
      }
      
      if (plot_pdf) {
        if(tag == "all") {
          fn = sprintf("%s/cumulative_time_dist_%d.pdf",fig_dir,i)
        } else {
          fn = sprintf("%s/cumulative_time_dist_%d_%s.pdf",fig_dir,i,tag)
        }
        pdf(fn,w = 5*fig_scale,h = 4*fig_scale,useDingbats = F)
        plot(c(xlim_min:xlim_max),host_dens[xlim_min:xlim_max],type = "l", main = paste0("Control ",i),lwd = lwd,
             xlab = "Transcriptional rank",ylab = "",ylim = c(0,1),col = "gray30",cex.main = 3,cex.lab = 2,cex.axis = 2)
        
        lines(c(xlim_min:xlim_max),control_dens[xlim_min:xlim_max],type = "l",lwd = lwd,col = "cornflowerblue")
        
        text(x = rep(x_pos_text,4), y = c(0.9,0.75,0.6,0.45),
             labels = c(sprintf("D = %1.1e",ks_control_host_chim$statistic),sprintf("p < %1.1e",ks_control_host_chim$p.value),
                        sprintf("N host = %d",N_host),sprintf("N Control = %d",N_control)),cex = 1.2)
        
        dev.off()
      } else {
        png(sprintf("%s/cumulative_time_dist_%s.png",fig_dir,chimera),w = 500*fig_scale,h = 400*fig_scale)
        plot(dev_time,host_dens,type = "l", main = chimera,lwd = lwd,
             xlab = "Developmental time",ylab = "",ylim = c(0,1),col = "gray30",cex.main = 3,cex.lab = 2,cex.axis = 2)
        
        lines(dev_time,control_dens,type = "l",lwd = lwd,col = "cornflowerblue")
        
        text(x = rep(x_pos_text,4), y = c(0.9,0.75,0.6,0.45),
             labels = c(sprintf("D = %1.1e",ks_control_host_chim$statistic),sprintf("p < %1.1e",ks_control_host_chim$p.value),
                        sprintf("N host = %d",N_host),sprintf("N control = %d",N_control)),cex = 1.2)
        dev.off()
      }
      
      
    }
 
  }
  
  return(list(d_ks = ks_statistic,p_value = p_value,delta_median = delta_median,chimera = chim_embryos))
}
