add_alpha = function(col, alpha) {
   return(rgb(t(col2rgb(col))/256, alpha=alpha))
}

persp_polygon = function(x, y, col, z, k_persp=0.5)
{
	polygon(y, 14-(x+z*k_persp), col=col, border=NA)
}

plot_persp_levels = function(smoo_y1, smoo_y2, col_persp, cols, lims, k_persp=0.5, level_step=20)
{	
	max_xi = length(smoo_y1[[1]])
	max_y = lims[2]
	min_y = lims[1]
	yrange = floor((max_y-min_y)*100)
	zs = matrix(ncol = yrange, nrow=max_xi,0)
	for(col in cols) {
		z = col_persp[col]
		y1 = smoo_y1[[col]]
		y2 = smoo_y2[[col]]
		if(sum(y1)>sum(y2)) {
			y2 = smoo_y1[[col]]
			y1 = smoo_y2[[col]]
		}

		idx = apply(cbind(1:max_xi,(y1-min_y)*100,(y2-min_y)*100), 1, function(y) {
								xi = y[1]
								y1bin = max(0,floor(y[2])-1)
								y2bin = floor(y[3])-1
								if((y2bin-y1bin)<2) {
									return(c()) 
								} else {
									return(xi + max_xi*((y1bin-2):(y2bin+2)))
								}
							})
		zs[unlist(idx)] = z
	}
	ys = seq(min_y, max_y, l=yrange)
	shades = colorRampPalette(c("cornsilk3", "darkgoldenrod2", "darkgoldenrod4", "black"))(200)
	for(xi in seq(1,max_xi,level_step)) {
		z = zs[xi,]
      z_lo = loess(z ~ ys, span=0.03)$fitted
		lev_cols = shades[pmax(pmin(floor(z_lo*100)+1,200),1)]
		sx = ys
		sy = 13-xi/100-z_lo*k_persp
		n = length(sx)
		segments(sx[-n],sy[-n],sx[-1],sy[-1], lty=1, lwd=0.3, col=lev_cols[-1])
	}

}

draw_sig_edge_3d = function(x1, x2, x2t, y1, y2, y2t, flow, col1, col2, col_alpha=0.8, clip_top, clip_bot, z1, z2) 
{
	x1 = x1
	y1 = y1
	dx = x2t - x1
	dy = y2t - y1

	y1t = y1+flow
	dxt = x2 - x1
	dyt = y2 - y1t

	dz = z2-z1

	col1 = col2rgb(col1)[,1]
	names(col1) = c("red","green","blue")
	col2 = col2rgb(col2)[,1]
	names(col2) = c("red","green","blue")
	res = 0.01
#	message("segs ", "x1 ", x1, " y1 ", y1, " dx ", dx, " dy ", dy, " y1t ", y1t, " dxt ", dxt, " dyt ", dyt)
	beta0 = plogis(0,loc=0.5,scale=0.2)
	beta_f = plogis(1,loc=0.5,scale=0.2)-plogis(0,loc=0.5, scale=0.2)
	for(r in seq(0,1,res)) {
			beta = (plogis(r,loc=0.5,scale=0.2)-beta0)/beta_f
			beta5 = (plogis(r+res,loc=0.5,scale=0.2)-beta0)/beta_f

			sx1 = x1+r*dx
			sy1 = y1+beta*dy
			sx2 = x1+(r+res)*dx
			sy2 = y1+beta5*dy

			sx1t = x1+r*dxt
			sy1t = y1t+beta*dyt
			sx2t = x1+(r+res)*dxt
			sy2t = y1t+beta5*dyt

			sz = z1+r*dz
			sz2 = z1+(r+res)*dz

			r_col = r

			sy1 = pmin(pmax(sy1, clip_bot[as.character(round(sx1,2))]),
										clip_top[as.character(round(sx1,2))])
			sy2 = pmin(pmax(sy2, clip_bot[as.character(round(sx2,2))]),
										clip_top[as.character(round(sx2,2))])
			sy1t = pmin(pmax(sy1t, clip_bot[as.character(round(sx1t,2))]),
										clip_top[as.character(round(sx1t,2))])
			sy2t = pmin(pmax(sy2t, clip_bot[as.character(round(sx2t,2))]),
										clip_top[as.character(round(sx2t,2))])
			rgb_r = col2["red"]*r_col+col1["red"]*(1-r_col)
			rgb_g = col2["green"]*r_col+col1["green"]*(1-r_col)
			rgb_b = col2["blue"]*r_col+col1["blue"]*(1-r_col)
			col = rgb(rgb_r/256, rgb_g/256, rgb_b/256, col_alpha)
			persp_polygon(c(sx1, sx2, sx2t,sx1t), 
								c(sy1, sy2, sy2t, sy1t), 
								col=col, 
								z = c(sz,sz2,sz2,sz))
	}
#		segments(x1,y1,x2,y2, 
}

plot_all_veins = function(ordered_cols,fig_dir="figs", plot_pdf = F,first_col="#635547", xlim=c(-4,12), col_persp=NULL, add_levels=T, fn="all") {
  
  mc = scdb_mc("sing_emb_wt10_recolored")
  
  mat = scdb_mat("sing_emb_wt10")
  
  md = mat@cell_metadata
  
  type_ag= table(mc@colors[mc@mc], md[names(mc@mc),"age_group"])
  
  type_agn = t(t(type_ag)/colSums(type_ag))
  
  mct = scdb_mctnetwork("sing_emb_wt10")
  type_flow = mctnetwork_get_type_flows(mct, 1,13)
  
  key = mc@color_key
  rownames(key) = key$group
  
  t1 = 1
  t2 = 12
  T_minflow_for_type = 0.005
  T_minflow = 6e-3
  k_space_z = 0.2
  
  cols = ordered_cols

  center_i = which(cols == first_col)

  top_front = rep(0,1+(t2-t1)*100)

  foc_agn = type_agn[cols,t1:t2]
  foc_agn[foc_agn < 1e-3] = 0

  smoo_y1 = list()
  smoo_y2 = list()
  x = t1:t2
  if(plot_pdf) {
      pdf(sprintf("%s/%s.pdf",fig_dir, fn),w=20,h=12,useDingbats = F)
  } else {
      png(sprintf("%s/%s.png",fig_dir, fn),w=3000,h=1600,bg = "transparent")
  }
    
  plot(0, ylim=c(t1 - 2,t2 + 2), xlim=xlim)
    
  y_expand = seq(1,2,l=1+(t2-t1)*100)**1
  ry_expand = rev(y_expand)
  prev_z = 0
  for(i in center_i:length(cols)) {
  		cur_col = cols[i]
      y = foc_agn[cur_col,t1:t2]
      ys = approx(x,y, seq(t1,t2,l=1+(t2-t1)*100))
      lo = loess(ys$y ~ ys$x,span=0.3)$fitted

      calpha = add_alpha(cur_col,0.8)
		if(i == center_i) {
      	smoo_y2[[cur_col]] = (lo+top_front)*y_expand
      	smoo_y1[[cur_col]] = (-lo+top_front)*y_expand
      	names(smoo_y1[[cur_col]]) = ys$x 
      	names(smoo_y2[[cur_col]]) = ys$x 
      	persp_polygon(c(ys$x,rev(ys$x)), 
						c((lo+top_front)*y_expand, (rev(top_front)+rev(-lo))*ry_expand), 
						col = calpha,
						z = rep(col_persp[cur_col],length(top_front)*2))
		} else {
      	smoo_y2[[cur_col]] = (2*lo+top_front)*y_expand
      	smoo_y1[[cur_col]] = top_front*y_expand
      	names(smoo_y1[[cur_col]]) = ys$x
      	names(smoo_y2[[cur_col]]) = ys$x
	      persp_polygon(c(ys$x,rev(ys$x)), 
					  c((2*lo+top_front)*y_expand, rev(top_front)*ry_expand), 
						col = calpha,
						z = rep(col_persp[cur_col], length(top_front)*2))
		}
		
		dz = max(0,col_persp[cur_col] - prev_z)
		prev_z = col_persp[cur_col]
		top_front = top_front + dz*k_space_z + loess(ifelse(lo>0, 0.05+0.5*max(lo)+2*lo, 0) ~ ys$x, span=0.7)$fitted
		if(i == center_i) {
			top_front = top_front / 2
			bot_front = -top_front
		}
  }

 prev_z = 0
 if(center_i != 1) {
  for(i in (center_i-1):1) {
  		cur_col = cols[i]
      y = foc_agn[cur_col,t1:t2]
      ys = approx(x,y, seq(t1,t2,l=1+(t2-t1)*100))
      lo = loess(ys$y ~ ys$x,span=0.3)$fitted

      calpha = add_alpha(cur_col,0.8)
     	smoo_y2[[cur_col]] = bot_front*y_expand
      smoo_y1[[cur_col]] = (-2*lo+bot_front)*y_expand
      names(smoo_y1[[cur_col]]) = ys$x
      names(smoo_y2[[cur_col]]) = ys$x
		#message("cur_col ", cur_col, " persp " , col_persp[cur_col])
      persp_polygon(c(ys$x,rev(ys$x)), 
				  c(bot_front*y_expand, (rev(bot_front)+rev(-2*lo))*ry_expand), 
					col = calpha,
					z = rep(col_persp[cur_col], length(top_front)*2))

		dz = max(0,col_persp[cur_col] - prev_z)
		prev_z = col_persp[cur_col]
		bot_front = bot_front - dz*k_space_z - loess(ifelse(lo>0, 0.05+0.5*max(lo)+2*lo, 0) ~ ys$x, span=0.7)$fitted
 	 }
  }
  if(add_levels) {
  	plot_persp_levels(smoo_y1, smoo_y2, col_persp, cols, xlim)
  }
 
  for(foc_i in 1:length(cols)) {
	foc_color = cols[foc_i]
   for(t in t1:(t2-1)) {
      flow = type_flow[[t]]
      max_i = length(cols)
      cum_y = smoo_y2[[foc_color]][as.character(t)]
		if(foc_i > 1) {
      for(i in 1:(foc_i-1)) {
        col_i = cols[i]
        fl = flow[foc_color,col_i]*2
        if(!is.null(fl) & length(fl) > 0 & fl > T_minflow) {
          calpha = add_alpha(col_i, 0.8)
          draw_sig_edge_3d(x1=t,
								x2 = t+1,
								x2t = t+1-2*fl-0.05,
                        y1 = smoo_y1[[foc_color]][as.character(t)], 
                        y2 = smoo_y2[[col_i]][as.character(t+1)], 
                        y2t =smoo_y2[[col_i]][as.character(round(t+1-2*fl-0.05,2))],
                        flow = fl,
                        col1=foc_color, col2=col_i,
								clip_top = smoo_y1[[foc_color]],
								clip_bot = smoo_y2[[col_i]],
								z1=col_persp[foc_color], z2=col_persp[col_i])
        }
      }
		}
		if(foc_i < max_i) {
      for(i in max_i:(foc_i+1)) {
        col_i = cols[i]
        fl = flow[foc_color,col_i]*2
        if(fl > T_minflow) {
          calpha = add_alpha(col_i, 0.8)
          draw_sig_edge_3d(x1=t,
								x2t = t+1,
								x2 = t+1-2*fl-0.05,
                        y1 = smoo_y2[[foc_color]][as.character(t)]-fl, 
                        y2t = smoo_y1[[col_i]][as.character(t+1)], 
                        y2 = smoo_y1[[col_i]][as.character(round(t+1-2*fl-0.05,2))],
                        flow = fl,
                        col1=foc_color, col2=col_i,
								clip_bot = smoo_y2[[foc_color]],
								clip_top = smoo_y1[[col_i]],
								z1=col_persp[foc_color], z2=col_persp[col_i])
        }
      }
		}
      #incoming
    }
	 }
    dev.off()
}


