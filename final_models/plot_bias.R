plot_bias <- function(stanfit) {
  
dem_color <- "#4589c3"
rep_color <- "#c34589"
landmark_color <- "#45c3be"

dark_gray_color <- "#222222"
theme <- function(...) ggplot2::theme(...)
transparent <- theme(panel.background = element_blank(), plot.background = element_blank())
axis_line_color <- dark_gray_color
axis_color <- theme(axis.line = element_line(color = axis_line_color))
axis_labs <- theme(axis.title = element_text(face = "bold", size = 13))
title_txt <- theme(plot.title = element_text(face = "bold", size = 14))
fat_axis <- theme(axis.line.x = element_line(size = 3, color = axis_line_color), 
                  axis.line.y = element_line(size = 0.5, color = axis_line_color))
h_lines <- theme(panel.grid.major = element_line(size = 0.10, linetype = 3, color = "turquoise4"),
                 panel.grid.major.x = element_blank())
v_lines <- theme(panel.grid.major = element_line(size = 0.25, linetype = 3, color = "turquoise4"),
                 panel.grid.major.y = element_blank())
no_yaxs <- theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())


plot_bias <- function(posterior, data, cred_lev = 0.95, legend = TRUE) {  
  prepare_plot_data <- function() {
    Mean <- apply(posterior, 2, mean)
    SD <- apply(posterior, 2, sd)
    Majority <- factor(data$majority, labels = c("Republican", "Democrat"))
    alpha <- 1 - (1 - cred_lev)/2
    z <- qnorm(alpha)
    out <- data.frame(Congress = data$congress, Majority, Bias = Mean, LB = Mean - z*SD, UB = Mean + z*SD)
    return(out)
  }
  assemble_graph <- function(df) {
    my_clrs <- scale_color_manual(values = c(rep_color, dem_color), name = "Majority Party")
    period_lab <- "Period hypothesized\n by Cox & Katz\n to show bias\n toward majority"
    my_lty <- scale_linetype_manual(values = 1, labels = period_lab, name = "")
    my_guides <- guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))
    x_scale <- scale_x_continuous(breaks = c(46, seq(50,100,10), 106), 
                                  labels = c("46th \n\n (1879-1881)", seq(50,100,10), "106th \n\n (1999-2001)"))
    rects1 <- annotate("rect", xmin=xmin, xmax=xmax, ymin=-0.175, ymax=-0.125, 
                       fill = c("black","gray25","black","gray25"))
    rects2 <- annotate("rect", xmin=xmin[c(2,4)], xmax=xmax[c(2,4)], ymin=-0.130, ymax=-0.125, fill = landmark_color)
    txts1 <- annotate("text", x = (xmin + xmax)/2, y = -0.1375,
                      label = landmark_names[,1], size = 3.5, color = "gray75")
    txts2 <- annotate("text", x = (xmin + xmax)/2, y = -0.16,
                      label = landmark_names[,2], size = 3.5, color = "gray75")
    g <- ggplot(df, aes(x = Congress, y = Bias, ymin = LB, ymax = UB))
    graph <- (g + my_lty + my_clrs + x_scale + ylab("Bias toward majority") +
                geom_segment(aes(linetype=""), x=xmin[2], xend=xmax[2], y=-.1725, yend=-0.1725, size=1, color=landmark_color) +
                rects1 + rects2 + txts1 + txts2 +
                geom_hline(yintercept = 0, color = "darkgray") +
                geom_linerange(aes(color = Majority), size = 1.75, alpha = 1) +
                geom_line(color = "black") +
                geom_point(color = "black", size = 1.5) +
                my_guides
    )
    if (!legend) graph <- graph + theme(legend.position = "none")
    return(graph)
  }
  df <- prepare_plot_data()
  graph <- assemble_graph(df)
  graph <- graph + theme_classic() %+replace% (fat_axis + axis_color)
  return(graph)
}

display_plots <- function(labs, graph1) {
  gr_ttl <- textGrob(labs["title"], gp = gpar(fontface="bold.italic", fontsize=17))
  gr_sub1 <- textGrob(labs$subs[1], gp = gpar(fontface="bold.italic", fontsize=14,lineheight=1))  
  gr_sub2 <- textGrob(labs$subs[2], gp = gpar(fontsize=13))
  space <- textGrob(" ")
  image <- arrangeGrob(gr_sub1,gr_sub2, graph1,
                       ncol=1, nrow=3, heights = c(1,.5,10))
  grid.arrange(image)
}



load("ck_data.RData")

library(rstan)
library(plyr)
library(ggplot2)
bias_posterior <- extract(stanfit, pars = "bias")[[1]]
congress_data <- ddply(ck_data, "congress", summarise, majority = mean(demmaj))


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
ttl <- c("\n Roll Calls in the US House of Representatives, 46th - 106th Congresses \n")
subs <- c("Estimated bias toward the majority party \n", "Posterior means and 95% credible intervals")


labs <- list(title = ttl, subs = subs)
graph1 <- plot_bias(bias_posterior, congress_data, cred_lev = 0.95, legend = T)
graph1

}



graph <- plot_bias(ck_cholesky)
ggsave(filename = "final_models/ck_bias.pdf", width = 9, height = 6)

library(shinyStan)
load("final_models/ck_cholesky.RData")
load("ck_data.RData")
y <- ck_data$majps # For pp checks
launch_shinystan(ck_cholesky)


library(gridExtra)
pdf(file = "ck_diagnostics.pdf", w = 9, h = 4)
grid.arrange(shinystan_n_eff, shinystan_mcse, shinystan_rhat, nrow = 1)
graphics.off()
save(shinystan_mcse, file = "final_models/ck_mcse_plot.RData", compress = "xz")
