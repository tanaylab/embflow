# fig 7 plots
source("scripts/generate_paper_figures/plot_3d_vein.r")

gen_fig_7_plots = function() {
  
  if(!dir.exists("figs/paper_figs/fig7")) {
    dir.create("figs/paper_figs/fig7")
  }
  
  fig_7a()
  fig_7b()
  fig_7c()

}


fig_7a = function(plot_pdf = T) {
  
  fig_dir = "figs/paper_figs/fig7"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  theiler_stage_order = c("7","8","9","10","10b","10c","11","11b","11c","11d","11e")
  
  mat = scdb_mat("sing_emb_wt10")
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  emb_age_time = unique(mat@cell_metadata[names(mc@mc),c("embryo","Theiler_Stage","developmental_time","age_group")])
  f = !is.na(emb_age_time$Theiler_Stage)
  emb_age_time = emb_age_time[f,]
  
  emb_age = split(emb_age_time$developmental_time,f = emb_age_time$Theiler_Stage)
  
  
  th_stage_median = round(tapply(emb_age_time$developmental_time,emb_age_time$Theiler_Stage,median),2)
  
  n_1 = (length(th_stage_median) + 1) %/% 2
  n_2 = length(th_stage_median) %/% 2
  odd_n = 2*c(1:n_1) - 1
  even_n = 2*c(1:n_2)
  
  if(plot_pdf) {
    pdf(sprintf("%s/fig_7a.pdf",fig_dir),w = 7,h = 5)
  } else {
    png(sprintf("%s/fig_7a.png",fig_dir),w = 1000,h = 500)
  }
  
  out = boxplot(emb_age[theiler_stage_order],horizontal = T,pch = 19,las = 2,xaxt = 'n',ylim = c(6.4,8.1))
  #axis(side = 1,at = c(6.5,6.75,7,7.25,7.5,7.75,8),labels = c(6.5,6.75,7,7.25,7.5,7.75,8)
  axis(side = 1,at = th_stage_median[odd_n],labels = th_stage_median[odd_n],cex.axis = 0.5)
  axis(side = 1,at = th_stage_median[even_n],labels = th_stage_median[even_n],cex.axis = 0.5)
  dev.off()
}


fig_7b = function(plot_pdf = T) {
  
  fig_dir = "figs/paper_figs/fig7/fig_7b"
  
  w = 300
  h = 600
  
  mat_id = "sing_emb_wt10"
  mc_id = "sing_emb_wt10_recolored"
  
  mat = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  
  col_to_rank = c(1:nrow(mc@color_key))
  names(col_to_rank) = mc@color_key$color
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  } 
  
  excluded_colors = c("#7F6874","#F6BFCB")
  #excluded_colors = c("#7F6874")
  
  
  cls = names(mc@mc)[!(mc@colors[mc@mc] %in% excluded_colors)]
  
  ct_vs_emb = table(mc@colors[mc@mc[cls]],mat@cell_metadata[cls,"transcriptional_rank"])
  ct_vs_emb = t(t(ct_vs_emb)/colSums(ct_vs_emb))
  
  col_ord = order(col_to_rank[rownames(ct_vs_emb)])
  
  ct_vs_emb = ct_vs_emb[col_ord,]
  
  emb_to_age_group = unique(mat@cell_metadata[names(mc@mc),c("age_group","transcriptional_rank")])
  emb_to_age_group = emb_to_age_group[order(emb_to_age_group$transcriptional_rank),]
  
  
  for (age in 1:max(emb_to_age_group$age_group)) {
    f = emb_to_age_group$transcriptional_rank[emb_to_age_group$age_group == age]
    if(plot_pdf) {
      pdf(sprintf("%s/age_%d.pdf",fig_dir,age),w= w/100,h= h/100)
    } else {
      png(sprintf("%s/age_%d.png",fig_dir,age),w= w,h= h)
    }
    par(mar = c(0.1,0.1,0.1,0.5))
    barplot(ct_vs_emb[,f], col = rownames(ct_vs_emb),las = 2,horiz = F,axes = F,axisnames = F)
    dev.off()
    
  }
  

  
  if(plot_pdf) {
    pdf(sprintf("%s/fig_7b.pdf",fig_dir),w= 13*w/100,h= h/100)
  } else {
    png(sprintf("%s/fig_7b.png",fig_dir),w= 13*w,h= h)
  }
  
  layout(matrix(c(1:13),nrow = 1,ncol = 13),w = rep(300,13),h = c(600))
  for (age in 1:max(emb_to_age_group$age_group)) {
    f = emb_to_age_group$transcriptional_rank[emb_to_age_group$age_group == age]
    par(mar = c(0.1,0.1,0.1,0.3))
    barplot(ct_vs_emb[,f], col = rownames(ct_vs_emb),las = 2,horiz = F,axes = F,axisnames = F)
  }
  dev.off()
  
  
  
} 




fig_7c = function() {
  
  # next follows the vein plot
  vein_par = read.table("data/external_data/fig7_vein_parameters.txt",sep = "\t",stringsAsFactors = F)
  ordered_cols = vein_par$color
  
  col_persp = vein_par$level
  names(col_persp) = vein_par$color
  
  # the following function is contained in scripts/generate_paper_figures/fig7_vein_all.r
  plot_all_veins(ordered_cols = ordered_cols,fig_dir = "figs/paper_figs/fig7",plot_pdf = T,col_persp = col_persp,fn = "fig_7c")
  
}

if(0) {
  size_infered_time = function() {
    
    fig_dir = "figs/paper_figs/fig7/size infered time"
    
    theiler_stage_order = c("7","8","9","10","10b","10c","11","11b","11c","11d","11e")
    
    mat = scdb_mat("sing_emb_wt10")
    mc = scdb_mc("sing_emb_wt10_recolored")
    
    emb_age_time = unique(mat@cell_metadata[names(mc@mc),c("embryo","Theiler_Stage","embryo_infered_time","age_group","morphology_rank","area","mouse_type")])
    
    f = (emb_age_time$mouse_type != "ICR") & !(is.na(emb_age_time$area))
    emb_age_time = emb_age_time[f,]
    f2 = emb_age_time$morphology_rank > 100 & log2(emb_age_time$area) < 16.3
    
    tmp = smooth.spline(x = emb_age_time[!f2,"morphology_rank"],log2(emb_age_time[!f2,"area"]),spar = 0.8)
    png(sprintf("%s/fit_morph_ranks_to_log_area.png",fig_dir))
    plot(emb_age_time[,"morphology_rank"],log2(emb_age_time[,"area"]),pch = 19,ylab = "log2(area)",xlab = "Morphology rank")
    lines(x = tmp$x,y = tmp$y)
    dev.off()
    
    emb_age_time2 = unique(mat@cell_metadata[names(mc@mc),c("embryo","Theiler_Stage","embryo_infered_time","age_group","morphology_rank","area","mouse_type")])
    f = !is.na(emb_age_time2$morphology_rank)
    emb_age_time2 = emb_age_time2[f,]
    emb_age_time2 = emb_age_time2[order(emb_age_time2$morphology_rank),]
    
    infered_time = (tmp$y - min(tmp$y))/(max(tmp$y) - min(tmp$y))*(8.1 - 6.4) + 6.4
    names(infered_time) = tmp$x
    
    morph_time = c(1:nrow(emb_age_time2))
    for(i in 1:max(tmp$x)) {
      if(i %in% tmp$x) {
        morph_time[i] = infered_time[as.character(i)]
      } else {
        last_ind = max(tmp$x[tmp$x < i])
        next_ind = min(tmp$x[tmp$x > i])
        morph_time[i] = ((i - last_ind)*infered_time[as.character(next_ind)] + (next_ind - i)*infered_time[as.character(last_ind)])/(next_ind - last_ind)
      }
    }
    
    for(i in (max(tmp$x) + 1):nrow(emb_age_time2)) {
      last_ind = max(tmp$x)
      scnd_last_ind = max(tmp$x[tmp$x < last_ind])
      
      morph_time[i] = morph_time[scnd_last_ind] + (morph_time[last_ind] - morph_time[scnd_last_ind])/(last_ind - scnd_last_ind)*(i - scnd_last_ind)
    }
    
    emb_age_time2$size_infered_time = morph_time
    
    tmp2 = smooth.spline(x = emb_age_time2[,"morphology_rank"],emb_age_time2[,"embryo_infered_time"],spar = 0.8)
    png(sprintf("%s/fit_morph_rank_to_transcriptional_infered_time.png",fig_dir))
    plot(emb_age_time2[,"morphology_rank"],emb_age_time2[,"embryo_infered_time"],pch = 19,ylab = "Transcriptional time",xlab = "Morphology rank")
    lines(x = tmp2$x,y = tmp2$y)
    dev.off()
    
    transcriptional_time = tmp$y
    
    png(sprintf("%s/compare_size_time_to_transcripional_time.png",fig_dir))
    plot(x = emb_age_time2$morphology_rank,y = tmp2$y,type = "l",ylab = "Infered time",xlab = "Morphology rank", col = "red")
    lines(x = emb_age_time2$morphology_rank,y = morph_time,col = "blue")
    legend(x = "topleft",legend = c("by size","by transcription"),col = c("blue","red"),lty = 1)
    dev.off()
    
    png(sprintf("%s/compare_size_time_to_transcripional_time2.png",fig_dir))
    plot(x = emb_age_time2$size_infered_time,y = emb_age_time2$embryo_infered_time,pch = 19,xlab = "Time infered from size", ylab = "Transcriptional time")
    dev.off()
    
    
    
    
    th_stage = split(emb_age_time2$size_infered_time,f = emb_age_time2$Theiler_Stage)
    th_stage_median = round(tapply(emb_age_time2$size_infered_time,emb_age_time2$Theiler_Stage,median),2)
    n_1 = (length(th_stage_median) + 1) %/% 2
    n_2 = length(th_stage_median) %/% 2
    odd_n = 2*c(1:n_1) - 1
    even_n = 2*c(1:n_2)
    
    png(sprintf("%s/size_infered_time.png",fig_dir),w = 1000,h = 500)
    boxplot(th_stage[theiler_stage_order],horizontal = T,pch = 19,las = 2,xaxt = 'n',ylim = c(6.4,8.1))
    axis(side = 1,at = th_stage_median[odd_n],labels = th_stage_median[odd_n],cex.axis = 0.8)
    axis(side = 1,at = th_stage_median[even_n],labels = th_stage_median[even_n],cex.axis = 0.8)
    dev.off()
    
    th_stage = split(tmp2$y,f = emb_age_time2$Theiler_Stage)
    th_stage_median = round(tapply(tmp2$y,emb_age_time2$Theiler_Stage,median),2)
    n_1 = (length(th_stage_median) + 1) %/% 2
    n_2 = length(th_stage_median) %/% 2
    odd_n = 2*c(1:n_1) - 1
    even_n = 2*c(1:n_2)
    
    png(sprintf("%s/transcription_infered_time.png",fig_dir),w = 1000,h = 500)
    boxplot(th_stage[theiler_stage_order],horizontal = T,pch = 19,las = 2,xaxt = 'n',ylim = c(6.4,8.1))
    axis(side = 1,at = th_stage_median[odd_n],labels = th_stage_median[odd_n],cex.axis = 0.8)
    axis(side = 1,at = th_stage_median[even_n],labels = th_stage_median[even_n],cex.axis = 0.8)
    dev.off()
    
    
    
  }
  
}
