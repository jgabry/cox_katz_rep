ck_pp_check <- function(stanfit, data, make_dir = TRUE) {
  require(rstan)
  require(arm)
  require(ggplot2)
  require(reshape2)
  
  
  graphics.off()
  print("Generating plots ...")
  
  stanfit_name <- deparse(substitute(stanfit))
  plot_dir <- paste0("plots/", stanfit_name)
  
  if (make_dir == TRUE) {
    print(paste("Creating directory", plot_dir))
    dir.create(plot_dir)
  }
  
  plot_dir <- paste0(plot_dir,"/")
  
  pdf(file = paste0(plot_dir, "ppcheck_plots.pdf"))
  
  y <- data$majps
  y_rep <- extract(stanfit, pars = "y_rep")[[1]]
  
  # residuals
  resids_rep <- y_rep
  for (i in 1:nrow(resids_rep)) {
    resids_rep[i,] <- y - resids_rep[i,]
  }
  
  # average y_rep
  avg_y_rep <- colMeans(y_rep)
  
  # average resids_rep
  avg_resids_rep <- colMeans(resids_rep)
  
  
  
  # histogram: residuals -------------------------------------------------
  hist(resids_rep,
       axes = FALSE, ylab = "", xlab = "",
       main = "Histogram of residuals")
  axis(1, lwd = 4)
  
  # binned residual plots for a sample of replicated data sets
  samp_size <- 4
  y_rep_samp <- y_rep[sample(nrow(y_rep), size = samp_size),]
  par(mfrow = c(2,2))
  for (i in 1:samp_size) {
    yr <- y_rep_samp[i, ]
    binnedplot(x = yr, y = y-yr, axes = FALSE,
               main = "", xlab = "Simulated values")
    axis(1, lwd = 4)
    axis(2, lwd = 4)
  }
  par(mfrow = c(1,1))
  mtext("Binned residual plots for four replicated data sets", side = 3, line = 2)
  
  
  
  # scatter: y vs avg_y_rep -------------------------------------------------
  plot(y, avg_y_rep, pch = 20, xlab = "Observed", ylab = "Avg. simulated", col = "purple", axes = FALSE)
  abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray35")
  axis(1, lwd = 4)
  axis(2, lwd = 4)
  
  
  
  # density: subset of simulated vs observed --------------------------------
  samp_size <- 20
  y_rep_samp <- y_rep[sample(nrow(y_rep), size = samp_size),]
  ymax <- max(c(density(y_rep_samp)$y, density(y)$y))
  
  plot(density(y_rep_samp[1,]), ylim = c(0, ymax), 
       axes = FALSE, xlab = "", ylab = "", 
       main = "Kernel density estimates (obs vs sample of reps)")
  for (i in 2:20) lines(density(y_rep_samp[i,]), col = "gray35")
  lines(density(y), col = "purple", lwd = 3)
  axis(1, lwd = 4)
  legend("topright", c("obs", "rep"), lwd = 2, col = c("purple", "black"), bty = "n")
  
  
  # scatter: avg_rep vs avg_resid_rep ---------------------------------------
  plot(avg_y_rep, avg_resids_rep, pch = 20, col = "purple",
       axes = FALSE, main = "Average replicated vs average residual",
       xlab = "Avg. replicated", ylab = "Avg. simulated residual")
  axis(1, lwd = 4)
  axis(2, lwd = 4)
  
  
  # Histogram: replications -------------------------------------------------
  hist(y_rep, freq = FALSE, 
       axes = FALSE, ylab = "", xlab = "",
       main = "Histogram of replications",
       # sub = list("Observed values plotted in purple along line y = 0 \n Dashed lines are medians", cex = 0.7),
       border = "white", col = "skyblue")
  # abline(v = median(y_rep), lty = 2, lwd = 1, col = "skyblue4")
  hy <- hist(y, plot = FALSE)
  lines(hy$mids, hy$density, col = "purple")
  # abline(v = median(y), lty = 2, col = "purple")
  axis(1, lwd = 4)
  legend("topright", c("obs", "rep"), lwd = 2, col = c("purple", "skyblue"), bty = "n")
  
  
  
  # histograms: test statistics ---------------------------------------------
  par(mfrow = c(2,2)) 
  funs <- c("mean", "sd", "min", "max")
  Tstats <- lapply(seq_along(funs), function(f) apply(y_rep, 1, funs[f]))
  names(Tstats) <- funs
  for (f in seq_along(funs)) {
    hist(Tstats[[f]], freq = FALSE, xlab = "", ylab = "", axes = FALSE, 
         border = "white", col = "skyblue", main = funs[f])
    abline(v = do.call(funs[f], args = list(y)), lwd = 3, col = "purple")
    axis(1, lwd = 4)
  }
  par(mfrow = c(1,1))
  mtext("Histograms of test statistics", side = 3, line = 2)
  
  
  
  # Time-series: data vs reps over time -------------------------------------
  
  # obs vs 20 replications
  period <- with(data, congress - 45)
  par(mfrow = c(2,1))
  sim_ids <- sample(1:nrow(y_rep), 20)
  y_rep_sample <- y_rep[sim_ids, ]
  
  plot(ts(tapply(y_rep_sample[1,], period, mean)), axes = FALSE,
       lwd = 0.75, col = "gray",
       main = "Observed data (purple) and \n 20 posterior predictive simulations (gray)", 
       xlab = "Period (Congress - 45)", ylab = "")
  for (i in 2:length(sim_ids)) {
    lines(ts(tapply(y_rep_sample[i,], period, mean)), 
          col = "gray", lwd = 0.75) 
  }
  lines(ts(tapply(y, period, mean)), col = "purple", lty = 2, lwd = 2)
  axis(1, lwd = 4)
  axis(2, lwd = 4)
  
  # plot obs (in purple) vs all sims
  plot(ts(tapply(y_rep[1,], period, mean)), axes = FALSE,
       lwd = 0.75, col = "gray",
       main = "Observed data (purple) and \n all posterior predictive simulations (gray)", 
       xlab = "Period (Congress - 45)", ylab = "")
  for (i in 2:nrow(y_rep)) {
    lines(ts(tapply(y_rep[i,], period, mean)), lwd = 0.75, col = "gray") 
  }
  lines(ts(tapply(y, period, mean)), col = "purple", lwd = 2, lty = 2)
  axis(1, lwd = 4)
  axis(2, lwd = 4)
  par(mfrow = c(1,1))
  
  
  
  # plots of estimated bias -------------------------------------------------  
  bias <- extract(stanfit, pars = "bias")[[1]]
  bias_m <- melt(bias, varnames = c("iter", "Congress"))
  
  fat_axes <- theme(axis.line = element_line(size = 3))
  my_theme <- theme_classic() %+replace% fat_axes

  
  majority <- with(data, tapply(demmaj, congress, mean))
  majority <- factor(majority, labels = c("Rep", "Dem"))
  bias_mean <- apply(bias,2,mean)
  bias_med <- apply(bias,2,median)
  bias_ub <- apply(bias,2,quantile,probs=0.975)
  bias_lb <- apply(bias,2,quantile,probs=0.025)
  inc_zero <- (bias_mean > 0 & bias_lb > 0) | (bias_mean < 0 & bias_ub < 0)
  inc_zero <- !inc_zero
  
   
  bias_dat <- data.frame(Congress = 1:length(bias_mean) + 45, bias_mean, bias_ub, bias_lb, majority, inc_zero)
  gg_bias_mean <- ggplot(bias_dat, aes(x = Congress, y = bias_mean, ymin = bias_lb, ymax = bias_ub)) + 
    geom_hline(yintercept = 0, color = "lightgray", lty = 2) + 
    geom_linerange(aes(color = majority), size = 3, alpha = 0.33) + 
    geom_line() +
    geom_point(size = 3) +
    scale_color_manual(values = c("red", "blue")) +
    labs(x = "Congress", y = "Bias") + 
    ggtitle("Posterior means and 95% intervals") +
    my_theme
  
  print(gg_bias_mean)
  
  graphics.off()
  
  
  
}

