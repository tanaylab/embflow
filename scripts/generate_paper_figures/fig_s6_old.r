library("Matrix")


plot_foxc1_foxc2_on_flow = function() {
  
  mct_id = "sing_emb_wt10"
  mat_id = "sing_emb_wt10"
  genes = c("Foxc1","Foxc2")
  mc = scdb_mc("sing_emb_wt10_recolored")
  mat = scdb_mat("sing_emb_wt10")
  mat_n = t(t(mat@mat[genes,names(mc@mc)])/colSums(mat@mat[,names(mc@mc)]))
  colors_ordered = mc@color_key$color
  
  cls_spl = split(names(mc@mc),mat@cell_metadata[names(mc@mc),"age_group"])
  
  reg = 3e-5
  
  legc_ls = lapply(genes,function(gene) {
    e_gc_mat = sapply(cls_spl,function(cls) {
      e_gc_tmp = tapply(mat_n[gene,cls],mc@mc[cls],mean)
      e_gc_t = rep(0,ncol(mc@e_gc))
      e_gc_t[as.numeric(names(e_gc_tmp))] = e_gc_tmp
      return(e_gc_t)
    })
    legc = log2(reg + e_gc_mat)
    return(legc)
  })
  
  
  mm_mctnetwork_plot_net(mct_id = mct_id,fn = "figs/paper_figs/fig_s6/foxc1_over_flows.png",colors_ordered = colors_ordered,mc_t_score = legc_ls[[1]],dy_ext = 0,dx_back = 0,w = 2200,h = 1200,tag = "Foxc1")
  mm_mctnetwork_plot_net(mct_id = mct_id,fn = "figs/paper_figs/fig_s6/foxc2_over_flows.png",colors_ordered = colors_ordered,mc_t_score = legc_ls[[2]],dy_ext = 0,dx_back = 0,w = 2200,h = 1200,tag = "Foxc2")
}

plot_color_scale_flow = function() {
  
  score_shades = colorRampPalette(c("lightgray", "gray", "darkgray", "lightpink", "pink", "red", "darkred"))(101)
  fig_fn = "figs/paper_figs/fig_s6/col_scale_flows_foxc1_foxc2.pdf"
  tgconfig::set_param(param = "mc_plot_device",value = "pdf",package = "metacell")
  plot_color_bar(vals = seq(0,1,length.out = 101),cols = score_shades,fig_fn = fig_fn,show_vals_ind = c(1,51,101))
  
  tgconfig::set_param(param = "mc_plot_device",value ="png",package = "metacell")
  

  
}


mm_mctnetwork_plot_net = function(mct_id, fn,
                                  mc_ord = NULL, colors_ordered = NULL,
                                  propogate=NULL,
                                  mc_t_score = NULL,
                                  edge_w_scale=5e-4,
                                  flow_thresh = 1e-4,
                                  w = 2000,h = 2000,
                                  mc_cex = 0.5,
                                  dx_back = 0.15, dy_ext = 0.4,
                                  sigmoid_edge = F, grad_col_edge = F,
                                  plot_mc_ids = F,miss_color_thresh = 0.5,
                                  func_deform=NULL,
                                  score_shades = colorRampPalette(c("lightgray", "gray", "darkgray", "lightpink", "pink", "red", "darkred"))(1000),
                                  tag = "")
{
  if(!is.null(propogate) | !is.null(mc_t_score)) {
    dx_back = 0
    dy_ext = 0
  }
  
  mct = scdb_mctnetwork(mct_id)
  if(is.null(mct)) {
    stop("cannot find mctnet object ", mct_id, " when trying to plot net flows")
  }
  net= mct@network
  mc = scdb_mc(mct@mc_id)
  if(is.null(mc)) {
    stop("cannot find mc object ", mct@mc_id, " matching the mc id in the mctnetwork object! db mismatch? recreate objects?")
  }
  if(is.null(mct@edge_flows) | sum(mct@edge_flows)==0) {
    stop("flows seems not to be initialized in mct id ", mct_id, " maybe rerun the mincost algorithm?")
  }
  names(mc@colors) = as.character(1:length(mc@colors))
  
  #color_ord = read.table("config/atlas_type_order.txt", h=T, sep="\t")
  #order MCs by type, mean age
  if(is.null(mc_ord)) {
    if(is.null(colors_ordered)) {
      stop("specify either mc_ord or color ord when plotting mctnet network")
    }
    mc_rank = mctnetwork_mc_rank_from_color_ord(mct_id, colors_ordered)
  } else {
    mc_rank = rep(-1,length(mc_ord))
    mc_rank[mc_ord] = c(1:length(mc_ord))
    names(mc_rank) = as.character(1:length(mc_rank))
  }
  
  mc_rank["-2"] = 0
  mc_rank["-1"] = length(mc_rank)/2
  
  #add growth mc
  
  f= net$flow > flow_thresh
  nn = net[f,]
  x1 = nn$time1
  x2 = nn$time2
  y1 = as.numeric(mc_rank[as.character(nn$mc1)])
  y2 = as.numeric(mc_rank[as.character(nn$mc2)])
  
  x1 = ifelse(nn$type1 == "growth", x1 + 0.3, x1)
  x2 = ifelse(nn$type2 == "growth", x2 + 0.3, x2)
  x1 = ifelse(nn$type1 == "norm_b" | nn$type1 == "extend_b",x1-dx_back,x1)
  x2 = ifelse(nn$type2 == "norm_b" | nn$type2 == "extend_b",x2-dx_back,x2)
  y1 = ifelse(nn$type1 == "src", max(y1)/2, y1)
  y2 = ifelse(nn$type2 == "sink", max(y2)/2, y2)
  y2 = ifelse(nn$type2 == "sink", NA, y2)
  y1 = ifelse(nn$type1 == "growth", y2+2.5, y1)
  y2 = ifelse(nn$type2 == "growth", y1+2.5, y2)
  y1 = ifelse(nn$type1 == "extend_b" | nn$type1 == "extend_f",y1+dy_ext, y1)
  y2 = ifelse(nn$type2 == "extend_f" | nn$type2 == "extend_b",y2+dy_ext, y2)
  
  if(!is.null(func_deform)) {
    x1 = ifelse(nn$type1 == "growth", NA, x1)
    x2 = ifelse(nn$type2 == "growth", NA, x2)
    y1 = ifelse(nn$type1 == "growth", NA, y1)
    y2 = ifelse(nn$type2 == "growth", NA, y2)
    min_x = min(c(x1,x2),na.rm=T)
    min_y = min(c(y1,y2),na.rm=T)
    range_x = (max(c(x1,x2),na.rm=T)-min_x)
    range_y = (max(c(y1,y2),na.rm=T)-min_y)
    ax = c((x1-min_x)/range_x, (x2-min_x)/range_x)
    ay = c((y1-min_y)/range_y, (y2-min_y)/range_y)
    atag = func_deform(ax, ay)
    n1 = length(x1)
    n = length(ax)
    x1 = atag[[1]][1:n1]
    x2 = atag[[1]][(n1+1):n]
    y1 = atag[[2]][1:n1]
    y2 = atag[[2]][(n1+1):n]
    
    if(0) {
      xy1 = func_deform((x1-min_x)/range_x, (y1-min_y)/range_y)
      xy2 = func_deform((x2-min_x)/range_x, (y2-min_y)/range_y)
      x1 = xy1[[1]]
      x2 = xy2[[1]]
      y1 = xy1[[2]]
      y2 = xy2[[2]]
    }
  }
  
  nn$mc1 = ifelse(nn$type1 == "src", nn$mc2, nn$mc1)
  
  png(fn, width = w,height = h)
  f_overflow = nn$type2=="extend_f" & nn$cost > 100
  f_underflow = nn$type2 == "norm_f" & nn$cost < -100
  #  nn$flow/(1e-8 + nn$capacity) < miss_color_thresh
  if(is.null(propogate) & is.null(mc_t_score)) {
    plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
    mc_rgb = col2rgb(mc@colors)/256
    f = nn$mc1>0 & nn$mc2 > 0
    m1 = as.numeric(nn$mc1[f])
    m2 = as.numeric(nn$mc2[f])
    seg_df = data.frame(x1 = x1[f], y1=y1[f], dx=x2[f]-x1[f], dy=y2[f]-y1[f], 
                        r1 = mc_rgb["red",m1],
                        r2 = mc_rgb["red",m2],
                        g1 = mc_rgb["green",m1],
                        g2 = mc_rgb["green",m2],
                        b1 = mc_rgb["blue",m1],
                        b2 = mc_rgb["blue",m2])
    
    for(alpha in seq(0,0.98,0.02)) {
      beta = alpha
      beta5 = alpha+0.02
      if(sigmoid_edge) {
        beta = plogis(alpha,loc=0.5,scale=0.1)
        beta5 = plogis(alpha+0.02,loc=0.5,scale=0.1)
      }
      sx1 = seg_df$x1+alpha*seg_df$dx
      sx2 = seg_df$x1+(alpha+0.02)*seg_df$dx
      sy1 = seg_df$y1+beta*seg_df$dy
      sy2 = seg_df$y1+beta5*seg_df$dy
      alpha_col = ifelse(grad_col_edge, alpha,0)
      rgb_r = seg_df$r2*alpha_col+seg_df$r1*(1-alpha_col)
      rgb_g = seg_df$g2*alpha_col+seg_df$g1*(1-alpha_col)
      rgb_b = seg_df$b2*alpha_col+seg_df$b1*(1-alpha_col)
      cols = rgb(rgb_r, rgb_g, rgb_b)
      segments(sx1, sy1, sx2, sy2, 
               col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
               lwd=pmin(nn$flow/edge_w_scale, 10))
    }
    #		segments(x1,y1,x2,y2, 
    #				col=ifelse(nn$type2=="growth", "black", mc@colors[nn$mc1]), 
    #				lwd=pmin(nn$flow/edge_w_scale, 10))
    f = f_overflow; segments(x1[f],y1[f],x2[f],y2[f], col="red", 
                             lwd=pmin(nn$flow[f]/edge_w_scale, 10))
    f = f_underflow; segments(x1[f],y1[f],x2[f],y2[f], col="blue", 
                              lwd=pmin((nn$capacity[f] - nn$flow[f])/edge_w_scale,10))
    #		points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=1)
  } else if(!is.null(mc_t_score)) {
    plot(c(x1,x2), c(y1,y2), pch=19, col="lightgray",cex=mc_cex,axes = F,xlab = "",ylab = "",main = tag)
    mc_t_score = pmax(mc_t_score, quantile(mc_t_score,0.03))
    mc_t_score = pmin(mc_t_score, quantile(mc_t_score,0.97))
    mc_t_score = mc_t_score-min(mc_t_score)
    mc_t_score = mc_t_score/max(mc_t_score)
    mc_t_score = floor(1+999*mc_t_score)
    f = nn$mc1>0 & nn$mc2>0 & nn$time1>0
    max_mc = nrow(mc_t_score)
    m1 = as.numeric(nn$mc1[f]) 
    m2 = as.numeric(nn$mc2[f]) 
    score1 = rep(1, nrow(nn))
    score2 = rep(1, nrow(nn))
    score1[f] = mc_t_score[m1+(nn[f,"time1"]-1)*max_mc]
    score2[f] = mc_t_score[m2+(nn[f,"time2"]-1)*max_mc]
    
    seg_df = data.frame(x1 = x1, y1=y1, dx=x2-x1, dy=y2-y1, 
                        score1= score1, dscore=score2-score1)
    
    for(alpha in seq(0,0.95,0.05)) {
      x1 = seg_df$x1+alpha*seg_df$dx
      x2 = seg_df$x1+(alpha+0.05)*seg_df$dx
      y1 = seg_df$y1+alpha*seg_df$dy
      y2 = seg_df$y1+(alpha+0.05)*seg_df$dy
      cols = score_shades[floor(seg_df$score1 + (alpha+0.025)*seg_df$dscore)]
      segments(x1, y1, x2, y2, 
               col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
               lwd=pmin(nn$flow/edge_w_scale, 10))
    }
    rect(seg_df$x1-0.12,seg_df$y1-0.5, seg_df$x1+0.12, seg_df$y1+0.5,
         col = mc@colors[nn$mc1], border=NA)

    #rect(seg_df$x1-0.12,seg_df$y1-0.5, seg_df$x1+0.12, seg_df$y1+0.5,
    #     col = alpha(mc@colors[nn$mc1],alpha =score1/1000), border=NA)
    #rect(seg_df$x2-0.12,seg_df$y2-0.5, seg_df$x2+0.12, seg_df$y2+0.5,
    #     col = mc@colors[nn$mc2], border=NA)
  } else {
    plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
    max_time = length(propogate)
    m1 = as.numeric(nn$mc1) 
    m2 = as.numeric(nn$mc2) 
    max_m = ncol(propogate[[1]])
    prop_flow = rep(0, nrow(nn))
    for(t in 1:max_time) {
      f = (nn$time1 == t) & nn$mc1>0 & nn$mc2>0
      prop_flow[f] = propogate[[t]][m1[f]+max_m*(m2[f]-1)]
    }
    segments(x1,y1,x2,y2, 
             col=ifelse(nn$type2=="growth", "black", mc@colors[nn$mc1]), 
             lwd=pmin(prop_flow/edge_w_scale, 10))
    points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=m_cex)
  }
  
  if(plot_mc_ids) {
    f1 = nn$type1!="growth" 
    text(x1[f1]-0.2,y1[f1], labels = nn$mc1[f1], cex=1)
    #	  text(c(x1[f1],x2[f2]),c(y1[f1],y2[f2]),labels = c(nn$mc1[f1],nn$mc2[f2]), cex=1)
  }
  
  dev.off()
}



gen_foxc1_foxc2_plus_one_tf_prediction = function() {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  col_to_ct = mc@color_key$group
  names(col_to_ct) = mc@color_key$color
  
  included_cell_types = c("Early nascent mesoderm","Late nascent mesoderm","Caudal mesoderm",
                          "Paraxial mesoderm","Cardiac crescent","Rostral mesoderm","Lateral & intermediate mesoderm","ExE mesoderm",
                          "Allantois","Haematoendothelial progenitors","Amnion/Chorion")
  
  mc_f = which(col_to_ct[mc@colors] %in% included_cell_types)
  
  legc = log2(mc@e_gc[,mc_f] + 3e-5)
  egc = mc@e_gc[,mc_f]
  legc_n = log2(mc@e_gc[,mc_f] + 3e-5) - log2(3e-5)
  legc_all = log2(mc@e_gc + 3e-5) - log2(3e-5)
  min_legc = apply(legc,1,min)
  max_legc = apply(legc,1,max)
  diff_legc = max_legc - min_legc
  
  f1 = (max_legc - min_legc > 2) & max_legc > -12
  
  feat_genes = rownames(legc)[f1]
  bad_genes = read.table("data/paper_data/gmods_wt10/bad_genes.txt",sep = "\t",stringsAsFactors = F)$x
  feat_genes = setdiff(feat_genes, bad_genes)
  ignore_tfs = c("Ccnd1","Naca","Hmga2","Psip1","Gtf2f2","Sall4","Basp1","Smad1","Actn1","Smad6")
  
  tf_f = read.table("data/all_tfs_list.txt",sep = "\t",stringsAsFactors = F)
  tf_f = tf_f$x
  
  marker = setdiff(feat_genes,tf_f)
  
  tf_list = intersect(feat_genes,tf_f)
  tf_list = setdiff(tf_list,ignore_tfs)
  
  marker = marker[diff_legc[marker] > 3]
  
  tf_list_f = setdiff(tf_list,c("Foxc1","Foxc2"))
  
  best_fit_1_tf_plus_foxc12 = function(targ) {
    r_sq = 0
    best_fit = c(0,0,0)
    for(i in 1:length(tf_list_f)) {
      tf1 = tf_list[i]
      m = summary(lm(legc[targ,] ~ legc[tf1,] + legc["Foxc1",] + legc["Foxc2",]))
      if(m$r.squared > r_sq) {
        best_fit = c(targ,tf1, m$r.squared)
        r_sq = m$r.squared
      }
    }
    return(best_fit)
  }
  
  best_pred_1_tf = sapply(marker,FUN =  function(g) {
    print(g)
    return(best_fit_1_tf_plus_foxc12(g))
  })
  
  best_pred_1_tf = t(best_pred_1_tf)
  colnames(best_pred_1_tf) = c("targ","tf1","m_rsq")
  
  write.table(best_pred_1_tf,file = "data/paper_data/fig5/tf_prediction/best_pred_1_tf_plus_foxc12.txt",sep = "\t")
  
  
}


plot_foxc12_plus_one_tf_vs_two_tf_pred = function() {
  
  fig_dir = "figs/paper_figs/fig_s6"
  
  best_pred = read.table("data/paper_data/fig5/tf_prediction/best_pred_2tfs.txt",sep = "\t",stringsAsFactors = F)
  best_pred_1_tf = read.table("data/paper_data/fig5/tf_prediction/best_pred_1_tf.txt",sep = "\t",stringsAsFactors = F)
  f1 = best_pred$m_rsq - best_pred_1_tf$m_rsq > 0.1
  
  best_pred_foxc12 = read.table("data/paper_data/fig5/tf_prediction/best_pred_1_tf_plus_foxc12.txt",sep = "\t",stringsAsFactors = F)
  rownames(best_pred_foxc12) = best_pred_foxc12$targ
  rownames(best_pred) = best_pred$targ
  
  plot(x = best_pred$m_rsq[f1],y = best_pred_foxc12$m_rsq[f1],pch = 19)
  f2 = best_pred_foxc12$m_rsq > 0.75 & best_pred$m_rsq - best_pred_foxc12$m_rsq <0.05
  text(x = best_pred$m_rsq[f2 &f1],y = best_pred_foxc12$m_rsq[f1 & f2],labels = best_pred$targ[f1 & f2])
  abline(a = 0,b = 1,lty = "dashed")
  
  
  
  
  
}