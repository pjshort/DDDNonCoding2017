# wrappers for commonly used visualization tools

sim_hist <- function(counts, observed, xlab = "Simulated Outcomes", ylab = "Simulation Frequency", main = "Comparing Observation to Simulation", col = "cyan", bin_width = 1){

  # makes a histogram of simulation outcomes and plots observed counts as dotted black line

  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  left_buffer = min(min(counts), observed)
  breaks = seq(-0.5, round(max(max(counts), observed)*1.2) + 0.5 + left_buffer/4, bin_width)
  h = hist(counts, xlab=xlab, ylab=ylab, main=main, breaks = breaks, col=col, xaxt="n")
  axis(side=1, at = breaks+0.5, labels=breaks + 0.5)
  abline(v=observed, col="black", lty=3, lw=5)
  
  empirical_p = sum(observed <= counts)/length(counts)
  print(sprintf("Empirical p-value for %s: %f", main, empirical_p))
}

density_hist_plot <- function(df, xfactor, fillfactor, title, xlabel, legend_title = "Data Set", binwidth=1){
  # xfactor and fillfactor should be factors of df
  ggplot(df, aes_string(x = xfactor, fill = fillfactor)) +
    geom_histogram(binwidth = 1, position = "dodge", aes(y = ..density..)) +
    geom_density(alpha = 0.50, kernel = "gaussian") +
    ggtitle(title) +
    xlab(xlabel) +
    ylab("Density") +
    guides(fill = guide_legend(override.aes = list(colour = NULL))) +
    scale_fill_discrete(name = legend_title) +
    theme(plot.title = element_text(size = 16),
          axis.text = element_text(size = 14))
}

# plot the poisson histogram against observed
poisson_bar_ggplot <- function(p, obs, main, min = 0, max = 1000){
  p = p[(min+1):(max+1)]
  counts = unlist(mapply(function(count, times) rep(count, times), seq(min, max), round(1000000*p)))
  sim_hist(counts, obs, main = main)
}

poisson_plot <- function(p, obs, main = "", min = 0, max = 1000, label_offset = 1.1) {
  pval = sum(p[seq(0, length(p) -1) > obs])
  p = p[(min+1):(max+1)]
  plot(seq(min, max), p, xlab = "Count", ylab = "Probability", main = main)
  abline(v=obs, col="black", lty=3, lw=5)
  text(x = obs*label_offset, y = max(p)*0.8, label = sprintf("p = %f\nobserved = %i\nexpected = %i", pval, obs, which.max(p) -1))
}
