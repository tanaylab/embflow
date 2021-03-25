library("metacell")
scdb_init("scrna_db/",force_reinit = T)

gen_fig_s3_plots = function() {
  
  if(!dir.exists("figs/paper_figs/fig_s3")) {
    dir.create("figs/paper_figs/fig_s3")
  }
  
  # both fig_s3a() and fig_s3b() require the file "data/fig_s3/mc_cluster_order.txt"
  # It contains a clustering of metacells based on the flows
  # Clustering of metacells is based on the network flow model
  # Users who want to redo this analysis should run the function cluster_metacells_by_flow(mct_id,K_cluster) 
  # and save the output in data/fig_s3/mc_cluster_order.txt
  fig_s3a()
  fig_s3a_color_scale()
  fig_s3b(T)
  
}

fig_s3b = function(plot_pdf = T) {
  fig_dir = "figs/paper_figs/fig_s3/fig_s3b"
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  mc = scdb_mc("sing_emb_wt10_recolored")

  mct_id = "sing_emb_wt10"
  K = 65
  marks = c("Foxc2","Lefty2","Tcf15","Nkx2-5","Myl4","Tagln","Pim2","Tcfap2a","Sox1","Grsf1","Dppa3","Etv2","Tal1","Cited4","Lefty1",
            "Noto","Foxa1","Ttr")
  
  mc_clust = read.table("data/fig_s3/mc_cluster_order.txt",sep = "\t",stringsAsFactors = F)

  plot_marks_along_clust(genes = marks,mc_clust = mc_clust,mc = mc,fig_dir = fig_dir,plot_pdf = plot_pdf)
}


fig_s3a = function(plot_pdf = T) {
  
  fig_dir = "figs/paper_figs/fig_s3"
  text_cex=1
  
  mct_id= "sing_emb_wt10"
  mct = scdb_mctnetwork(mct_id)
  
  cls_ord = read.table("data/fig_s3/mc_cluster_order.txt",sep = "\t",stringsAsFactors = F)
  mc_ord = cls_ord$mc
  
  K = length(unique(cls_ord$cluster))
  
  clst_flows = mctnetwork_clust_flows(mct_id = mct_id,K = K)
  
  cmat = clst_flows$cmat
  fclst = cls_ord$cluster[order(cls_ord$mc)]
  
  mc = scdb_mc(mct@mc_id)
  
  
  shades = colorRampPalette(c("darkblue", "blue","white", "red", "yellow"))(1000)
  if(plot_pdf) {
    pdf(sprintf("%s/fig_s3a.pdf",fig_dir), w=16, h=16)
    layout(matrix(c(1,2),nrow=2),heights=c(nrow(cmat)*0.09+0.50, 3))
    fig_scl = 1
  } else {
    fig_scl = 2
    png(sprintf("%s/fig_s3a.png",fig_dir), w=1600*fig_scl, h=1600*fig_scl)
    layout(matrix(c(1,2),nrow=2),heights=c(nrow(cmat)*9+50, 300)*fig_scl)
    
    
  }
  
  n_mc = nrow(cmat)
  
  
  par(mar=c(0,5,4,5))
  image(cmat[mc_ord, mc_ord], zlim=c(-1,1), col=shades, yaxt='n', xaxt='n')
  N = length(mc_ord)
  
  mc_x = 1:length(mc_ord)
  names(mc_x) = 1:N
  mc_x[mc_ord] = 1:N
  for(i in 1:K) {
    abline(h=max(-0.5+mc_x[fclst == i])/(N-1),lwd = fig_scl)
    abline(v=max(-0.5+mc_x[fclst == i])/(N-1),lwd = fig_scl)
  }
  cl_x = tapply(mc_x, fclst, mean)
  cl_max = tapply(mc_x, fclst, max)
  
  mtext(1:K, side = 3, at=cl_x/N, las=1, cex=fig_scl)
  par(mar=c(3,5,0,5))
  image(as.matrix(mc_ord,nrow=1), col=mc@colors, yaxt='n', xaxt='n')
  dev.off()
  
}

fig_s3a_color_scale = function() {
  
  # plot color scale
  shades = colorRampPalette(c("darkblue", "blue","white", "red", "yellow"))(101)
  
  vals = seq(-1,1,length.out = 101)
  cols = shades
  show_vals_ind =  c(1,51,101)
  
  
  pdf(file = "figs/paper_figs/fig_s3/fig_s3a_color_scale.pdf",useDingbats = F)
  plot.new()
  plot.window(xlim=c(0,100), ylim=c(0, length(cols) + 3))
  rect(7, 1:length(cols), 17, 1:length(cols) + 1, border=NA, col=cols)
  rect(7, 1, 17, length(cols)+1, col=NA, border = 'black')
  
  text(19, (1:length(cols))[show_vals_ind] + 0.5, labels=vals[show_vals_ind], pos=4)
  dev.off()
  
  
  
}



plot_marks_along_clust = function(genes,mc_clust,mc,fig_dir,plot_pdf = F,additional_horizontal_line = NULL) {
  
  fig_scl = 3
  
  mc_clust = mc_clust[order(mc_clust$mc_rank),]
  
  n_clust = max(mc_clust$cluster)
  
  mc_ord = mc_clust$mc
  
  fclst = mc_clust$cluster[order(mc_clust$mc)]
  
  legc = log2(mc@e_gc+1e-5)
  
  mc_x = 1:length(mc_ord)
  names(mc_x) = 1:length(mc_ord)
  mc_x[mc_ord] = 1:length(mc_ord)
  
  cl_x = tapply(mc_x, fclst, mean)
  cl_max = tapply(mc_x, fclst, max)
  
  
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  for(g in genes) {
    
    if(plot_pdf) {
      pdf(sprintf("%s/%s.pdf", fig_dir, g),w=10,h=3,useDingbats = F)
      #svg(sprintf("%s/%s.svg", fig_dir, g),w=10,h=3)
      fig_scl = 1
    } else {
      png(sprintf("%s/%s.png", fig_dir, g),w=1000*fig_scl,h=300*fig_scl)
    }
    plot(1:length(mc_ord), legc[g, mc_ord], pch=19, col=mc@colors[mc_ord], ylab="",xaxt='n',xlab = "",cex = fig_scl,cex.axis = fig_scl)
    mtext(1:n_clust, at=cl_x,side=1, las=2,cex = 0.6*fig_scl)
    
    if(!is.null(additional_horizontal_line)) {
      abline(h = additional_horizontal_line,lwd = fig_scl*3,lty = "dashed")
    }
    abline(v=cl_max+0.5,lwd = fig_scl)
    
    grid(lwd = fig_scl)
    dev.off()
  }
  
}







