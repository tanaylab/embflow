library("zoo")


chim_plot_ct_frequency_per_ct = function(plot_pdf = F) {
  
  mat_nm1 = "foxc_chim_wt10"
  mat_nm2 = "control_chim_wt10"
  
  rank_to_time = read.table(file = "data/revision/wt10_transcriptional_rank_developmental_time.txt",stringsAsFactors = F,h = T,sep = "\t")
  dev_time = rank_to_time$developmental_time
  dev_time = c(1:length(dev_time))
  
  age_field_host = "best_rank_host"
  age_field_ko = "best_rank_ko"
  age_field_control = "best_rank_control"
  
  
  fig_dir = sprintf("%s/fig6/ct_freq",.scfigs_base)
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  roll_width = 4
  
  load(file = sprintf("data/revision/%s/color_annotation/cmp_annot.Rda",mat_nm1))
  cmp_annot1 = cmp_annot
  load(file = sprintf("data/revision/%s/color_annotation/cmp_annot.Rda",mat_nm2))
  cmp_annot2 = cmp_annot
  
  query_cls_col1 = cmp_annot1$query_cls_col
  query_cls1 = names(query_cls_col1)
  query_cls_col2 = cmp_annot2$query_cls_col
  query_cls2 = names(query_cls_col2)
  
  
  mat1 =  scdb_mat(mat_nm1)
  mat2 =  scdb_mat(mat_nm2)
  mc_wt = scdb_mc("sing_emb_wt10_recolored")
  cgraph = scdb_cgraph(mat_nm1)
  col_to_ct = mc_wt@color_key$group
  names(col_to_ct) = mc_wt@color_key$color
  excluded_colors = c("#F6BFCB","#7F6874")
  included_colors = setdiff(unique(mc_wt@colors),excluded_colors)
  
  df_chim1 = read.table(sprintf("data/revision/%s/time_match/time_match_summary.txt",mat_nm1),sep = "\t",h = T,stringsAsFactors = F)
  rownames(df_chim1) = df_chim1$embryo
  df_chim2 = read.table(sprintf("data/revision/%s/time_match/time_match_summary.txt",mat_nm2),sep = "\t",h = T,stringsAsFactors = F)
  rownames(df_chim2) = df_chim2$embryo
  
  f = !(df_chim1$embryo %in% c("6D_m2e3","6D_m3e1","6D_m3e7"))
  chim_embryos1 = df_chim1$embryo[f]
  chim_embryos1 = chim_embryos1[order(df_chim1[chim_embryos1,"best_rank_host"])]
  chim_embryos2 = df_chim2$embryo[order(df_chim2[,"best_rank_host"])]
  
  ko_cls = query_cls1[( mat1@cell_metadata[query_cls1,"cell_type"] ==  "KO" ) & ( query_cls_col1[query_cls1] %in% included_colors )]
  control_cls = query_cls2[( mat2@cell_metadata[query_cls2,"cell_type"] ==  "control" ) & ( query_cls_col2[query_cls2] %in% included_colors )]
  host_cls1 = query_cls1[( mat1@cell_metadata[query_cls1,"cell_type"] %in% c("host") ) & ( query_cls_col1[query_cls1] %in% included_colors )]
  host_cls2 = query_cls2[( mat1@cell_metadata[query_cls2,"cell_type"] %in% c("host") ) & ( query_cls_col2[query_cls2] %in% included_colors )]
  wt10_cls = names(mc_wt@mc)[mc_wt@colors[mc_wt@mc] %in% included_colors]
  
  emb_wt_age = unique(mat1@cell_metadata[wt10_cls,c("transcriptional_rank","age_group")])
  emb_wt_age = emb_wt_age[order(emb_wt_age$transcriptional_rank),]
  
  wt10_emb_vs_ct = table(mat1@cell_metadata[wt10_cls,"transcriptional_rank"],mc_wt@colors[mc_wt@mc[wt10_cls]])
  ko_emb_vs_ct = table(mat1@cell_metadata[ko_cls,"embryo"],query_cls_col1[ko_cls])
  control_emb_vs_ct = table(mat2@cell_metadata[control_cls,"embryo"],query_cls_col2[control_cls])
  host_emb_vs_ct1 = table(mat1@cell_metadata[host_cls1,"embryo"],query_cls_col1[host_cls1])
  host_emb_vs_ct2 = table(mat2@cell_metadata[host_cls2,"embryo"],query_cls_col2[host_cls2])
  ko_emb_vs_ct = ko_emb_vs_ct[chim_embryos1,]
  control_emb_vs_ct = control_emb_vs_ct[chim_embryos2,]
  host_emb_vs_ct1 = host_emb_vs_ct1[chim_embryos1,]
  host_emb_vs_ct2 = host_emb_vs_ct2[chim_embryos2,]
  
  wt10_emb_vs_ct_n = wt10_emb_vs_ct/rowSums(wt10_emb_vs_ct[,included_colors])
  
  # next calculate moving  90% moving average window for every cell type
  
  mov_mean = apply(wt10_emb_vs_ct_n,2,function(freq_ct) {
    
    n = length(freq_ct)
    
    freq_ct = c(freq_ct,freq_ct[(n - roll_width + 1):n])
    
    freq_mean = rollmean(x = freq_ct,k = 2*roll_width+1)
    return(freq_mean)
  })
  
  mov_sd = apply(wt10_emb_vs_ct_n,2,function(freq_ct) {
    
    freq_ct = c(freq_ct,freq_ct[(length(freq_ct) - roll_width + 1):length(freq_ct)])
    
    freq_sd = rollapply(data = freq_ct,width = 2*roll_width+1,sd)
    
    return(freq_sd)
  })
  
  upper_sd = mov_mean + mov_sd
  lower_sd = mov_mean - mov_sd
  
  
  upper_lim = apply(wt10_emb_vs_ct_n,2,function(freq_ct) {
    
    freq_ct = c(freq_ct,freq_ct[(length(freq_ct) - roll_width + 1):length(freq_ct)])
    
    upper_freq = rollapply(data = freq_ct,width = 2*roll_width+1,function(v) {
      return(quantile(v,0.95))
    })
    
    return(upper_freq)
  })
  
  lower_lim = apply(wt10_emb_vs_ct_n,2,function(freq_ct) {
    
    freq_ct = c(freq_ct,freq_ct[(length(freq_ct) - roll_width + 1):length(freq_ct)])
    
    lower_freq = rollapply(data = freq_ct,width = 2*roll_width+1,function(v) {
      return(quantile(v,0.05))
    })
    
    return(lower_freq)
  })
  
  #rownames(upper_lim) = c((roll_width+1):(153 - roll_width))
  #rownames(lower_lim) = c((roll_width+1):(153 - roll_width))
  
  rownames(upper_lim) = c((roll_width+1):153)
  rownames(lower_lim) = c((roll_width+1):153)
  rownames(upper_sd) = c((roll_width+1):153)
  rownames(lower_sd) = c((roll_width+1):153)
  rownames(mov_mean) = c((roll_width+1):153)
  
  ko_emb_vs_ct_n = ko_emb_vs_ct/rowSums(ko_emb_vs_ct[,intersect(included_colors,colnames(ko_emb_vs_ct))])
  control_emb_vs_ct_n = control_emb_vs_ct/rowSums(control_emb_vs_ct[,intersect(included_colors,colnames(control_emb_vs_ct))])
  host_emb_vs_ct1_n = host_emb_vs_ct1/rowSums(host_emb_vs_ct1[,intersect(included_colors,colnames(host_emb_vs_ct1))])
  host_emb_vs_ct2_n = host_emb_vs_ct2/rowSums(host_emb_vs_ct2[,intersect(included_colors,colnames(host_emb_vs_ct2))])
  
  emb_vs_ct_ls = list(ko = ko_emb_vs_ct_n,host1 = host_emb_vs_ct1_n,wt10 = wt10_emb_vs_ct_n,host2 = host_emb_vs_ct2_n,control = control_emb_vs_ct_n)
  
  
  
  min_rank_plot = 105
  
  xlim_min = dev_time[min_rank_plot]
  #xlim_min = 7.6
  xlim_max = dev_time[153]
  
  for(cell_type in included_colors) {
    print(cell_type)
    
    
    ct_freq_ls = lapply(emb_vs_ct_ls,FUN = function(emb_vs_ct) {
      if(cell_type %in% colnames(emb_vs_ct)) {
        ct_freq = emb_vs_ct[,cell_type]
      } else {
        ct_freq = rep(0,nrow(emb_vs_ct))
        names(ct_freq) = rownames(emb_vs_ct)
      }
      return(ct_freq)
    })
    
    #calculate 90% moving average window
    wt10_freq = ct_freq_ls$wt10
    
    #x_ranks = c(min_rank_plot:(153-roll_width))
    x_ranks = c(min_rank_plot:(153))
    #upper_freq = upper_lim[as.character(c(min_rank_plot:(153-roll_width))),cell_type]
    #lower_freq = lower_lim[as.character(c(min_rank_plot:(153-roll_width))),cell_type]
    upper_freq = upper_sd[as.character(c(min_rank_plot:(153))),cell_type]
    lower_freq = lower_sd[as.character(c(min_rank_plot:(153))),cell_type]
    
    
    ylim_max = max(ct_freq_ls$ko,ct_freq_ls$host1,ct_freq_ls$host2,ct_freq_ls$control,ct_freq_ls$wt10)
    
    cell_type_nm = gsub("/","_",col_to_ct[cell_type])
    
    if(plot_pdf) {
      pdf(sprintf("%s/%s.pdf",fig_dir,cell_type_nm),w = 7,h = 5.6,useDingbats = F)
      par(mar = c(6,6,6,4))
      plot(x = dev_time[emb_wt_age$transcriptional_rank],y = ct_freq_ls$wt10,ylim = c(0,1.2*ylim_max),
           pch = 19,cex = 1,col = "gray80",main = col_to_ct[cell_type],ylab = "Fraction",xlab = "Transcriptional rank",cex.main = 3,
           cex.lab = 2,cex.axis = 2,xaxs = 'i',yaxs= 'i',xlim = c(xlim_min,xlim_max))
      
      
      polygon(x = c(dev_time[x_ranks],dev_time[rev(x_ranks)]),y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
      
      lines(x = dev_time[x_ranks],y = mov_mean[as.character(x_ranks),cell_type])
      
      points(x = dev_time[emb_wt_age$transcriptional_rank],y = ct_freq_ls$wt10,
             pch = 19,cex = 1,col = "gray50")
      
      segments(x0 = dev_time[df_chim1[names(ct_freq_ls$host1),age_field_host]],y0 = ct_freq_ls$host1,x1 = dev_time[df_chim1[names(ct_freq_ls$ko),age_field_host]],y1 = ct_freq_ls$ko,lwd = 1)
      
      segments(x0 = dev_time[df_chim2[names(ct_freq_ls$host2),age_field_host]],y0 = ct_freq_ls$host2,x1 = dev_time[df_chim2[names(ct_freq_ls$control),age_field_host]],y1 = ct_freq_ls$control,lwd = 1)
      
      
      points(x = dev_time[df_chim1[names(ct_freq_ls$host1),age_field_host]],y = ct_freq_ls$host1, pch = 19,cex = 3,col = "black")
      points(x = dev_time[df_chim1[names(ct_freq_ls$ko),age_field_host]],y = ct_freq_ls$ko, pch = 19,cex = 3,col = "#83c26d")
      
      points(x = dev_time[df_chim2[names(ct_freq_ls$host2),age_field_host]],y = ct_freq_ls$host2, pch = 17,cex = 3,col = "black")
      points(x = dev_time[df_chim2[names(ct_freq_ls$control),age_field_host]],y = ct_freq_ls$control, pch = 17,cex = 3,col = "cornflowerblue")
      
      dev.off()
      
    } else {
      
      png(sprintf("%s/%s.png",fig_dir,cell_type_nm),w = 800,h = 600)
      par(mar = c(6,6,6,4))
      plot(x = dev_time[emb_wt_age$transcriptional_rank],y = ct_freq_ls$wt10,ylim = c(0,1.2*ylim_max),
           pch = 19,cex = 1,col = "gray80",main = col_to_ct[cell_type],ylab = "Fraction",xlab = "Transcriptional rank",cex.main = 3,
           cex.lab = 3,xaxs = 'i',yaxs= 'i',xlim = c(xlim_min,xlim_max))
      
      
      polygon(x = c(dev_time[x_ranks],dev_time[rev(x_ranks)]),y = c(upper_freq,rev(lower_freq)),col = "gray80",border = NA)
      
      lines(x = dev_time[x_ranks],y = mov_mean[as.character(x_ranks),cell_type])
      
      points(x = dev_time[emb_wt_age$transcriptional_rank],y = ct_freq_ls$wt10,
             pch = 19,cex = 1,col = "black")
      
      segments(x0 = dev_time[df_chim1[names(ct_freq_ls$host1),age_field_host]],y0 = ct_freq_ls$host1,x1 = dev_time[df_chim1[names(ct_freq_ls$ko),age_field_host]],y1 = ct_freq_ls$ko,lwd = 1)
      
      segments(x0 = dev_time[df_chim2[names(ct_freq_ls$host2),age_field_host]],y0 = ct_freq_ls$host2,x1 = dev_time[df_chim2[names(ct_freq_ls$control),age_field_host]],y1 = ct_freq_ls$control,lwd = 1)
      
      
      points(x = dev_time[df_chim1[names(ct_freq_ls$host1),age_field_host]],y = ct_freq_ls$host1, pch = 19,cex = 3,col = "black")
      points(x = dev_time[df_chim1[names(ct_freq_ls$ko),age_field_host]],y = ct_freq_ls$ko, pch = 19,cex = 3,col = "#83c26d")
      
      points(x = dev_time[df_chim2[names(ct_freq_ls$host2),age_field_host]],y = ct_freq_ls$host2, pch = 17,cex = 3,col = "black")
      points(x = dev_time[df_chim2[names(ct_freq_ls$control),age_field_host]],y = ct_freq_ls$control, pch = 17,cex = 3,col = "cornflowerblue")
      
      legend(x = "topleft",legend = c("wt","host KO","host control","KO","control"),pch = c(19,19,17,19,17),col = c("black","black","black","#83c26d","cornflowerblue"))
      
      dev.off()
    }
    
    
  }
  
}


chim_plot_legend = function() {
  
  pdf("figs_revision/chimera_foxc_control/ct_freq/legend_cell_types.pdf",useDingbats = F)
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend(x = "center",legend = c("WT","Host KO","Host Control","KO","Control"),pch = c(19,19,17,19,17),col = c("gray50","black","black","#83c26d","cornflowerblue"))
  dev.off()
  
  
}