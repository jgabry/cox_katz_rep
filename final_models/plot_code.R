grid_arrange_shared_legend <- function(...) {
  # https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
  require(gridExtra)
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}



make_plot <- function() {
  congress <- data.frame(number = 46:106)
  congress <- mutate(congress, 
                     start = 1878 + seq(1,122,2), 
                     end = start + 2,
                     pre_Reed = (end <= 1890),
                     post_Czar = (end %in% 1890:1910),
                     post_Cannon = (end %in% 1912:1961),
                     post_Packing = (end %in% 1962:2001)
  )
  landmarks <- c(1890, 1910, 1961)
  ll <- length(landmarks)
  xmin <- congress$number[1]
  xmax <- congress$number[nrow(congress)]
  for(i in 1:ll){
    xmin <- c(xmin, congress$number[which.max(which(congress$end <= landmarks[i]))])
  }
  xmax <- c(xmin[-1], xmax)
  landmark_names <- cbind(c("pre-", rep("post-",3)), c("Reed", "Czar", "Cannon", "Packing"))
  
  
  
  dem_color <- "#4589c3"
  rep_color <- "#c34589"
  landmark_color <- "#45c3be"
  dark_gray_color <- "#222222"
  
  axis_line_color <- dark_gray_color
  axis_color <- theme(axis.line = element_line(color = axis_line_color))
  axis_labs <- theme(axis.title = element_text(face = "bold", size = 13))
  title_txt <- theme(plot.title = element_text(face = "bold", size = 14))
  fat_axis <- theme(axis.line.x = element_line(size = 3, color = axis_line_color), 
                    axis.line.y = element_line(size = 0.5, color = axis_line_color))
  lgnd_bot <- theme(legend.position = "bottom", legend.text = element_text(size = 12))
  
  
  my_clrs <- scale_color_manual(values = c(rep_color, dem_color), name = "Majority Party")
  period_lab <- "Period hypothesized by Cox & Katz to show bias toward majority"
  my_lty <- scale_linetype_manual(values = 1, labels = period_lab, name = "")
  my_guides <- guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))
  x_scale <- scale_x_continuous(breaks = c(46, seq(50,100,10), 106), 
                                labels = c("46th \n\n (1879-1881)", seq(50,100,10), "106th \n\n (1999-2001)"))
  y_scale <- scale_y_continuous(limits = c(-.2, .3))
  rects1 <- annotate("rect", xmin=xmin, xmax=xmax, ymin=-0.175, ymax=-0.125, 
                     fill = c("black","gray25","black","gray25"))
  rects2 <- annotate("rect", xmin=xmin[c(2,4)], xmax=xmax[c(2,4)], ymin=-0.130, ymax=-0.125, fill = landmark_color)
  txts1 <- annotate("text", x = (xmin + xmax)/2, y = -0.1375,
                    label = landmark_names[,1], size = 3.5, color = "gray75")
  txts2 <- annotate("text", x = (xmin + xmax)/2, y = -0.16,
                    label = landmark_names[,2], size = 3.5, color = "gray75")
  g <- ggplot(df, aes(x = Congress, y = Bias, ymin = LB, ymax = UB))
  graph <- (g + my_lty + my_clrs + x_scale + y_scale + ylab("Bias toward majority") +
              geom_segment(aes(linetype=""), x=xmin[2], xend=xmax[2], y=-.1725, yend=-0.1725, size=1, color=landmark_color) +
              rects1 + rects2 + txts1 + txts2 +
              geom_hline(yintercept = 0, color = "darkgray") +
              geom_linerange(aes(color = Majority), size = 1.75, alpha = 1) +
              geom_line(color = "black") +
              geom_point(color = "black", size = 1.5) +
              my_guides
  )
  
  graph <- graph + theme_classic() %+replace% (fat_axis + axis_color + lgnd_bot)
  
  graph
}
