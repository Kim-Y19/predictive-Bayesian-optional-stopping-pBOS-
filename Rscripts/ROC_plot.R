# Load necessary libraries
library(writexl)
library(ggplot2)
library(GGally)
library(dplyr)
library(pROC)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# Load the data
file_load_name <- "result_20241026_statistics_group400_Regressionyes_prior_var1_offsetmean_0_n0_10_v0_10_centermuTRUE_centervarTRUE_informativemuTRUE_informativevarTRUE_weaklyFALSE_statisticstandard_NEW1025__group_n_1"
load(paste0("Results/", file_load_name, ".RData"))

# Set color palette and line types
colors <- brewer.pal(8, "Dark2")

# Set theme for plots
ptsize <- 16
theme_set(theme_bw())
theme_update(
  axis.text = element_text(size = ptsize, colour = "black", family = "serif"),
  axis.line = element_line(colour = "black", size = 0.25),
  axis.ticks = element_line(colour = "black", size = 0.25),
  axis.title.y = element_text(margin = margin(t = 0, r = 0.0, b = 0, l = 0.2, unit = 'cm')),
  legend.key.width = unit(1.5, "cm"),
  legend.key.height = unit(0.4, "cm"),
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing = unit(0, "cm"),
  legend.position = "bottom",
  legend.text = element_text(size = ptsize, colour = "black", family = "serif"),
  legend.title = element_text(size = ptsize, colour = "black", family = "serif"),
  strip.background.x = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  plot.subtitle = element_text(size = ptsize, colour = "black", family = "serif", face = "plain", hjust = 0),
  plot.title = element_text(size = ptsize, colour = "black", family = "serif", face = "plain", hjust = 0.5, vjust = 0),
  text = element_text(size = ptsize, colour = "black", family = "serif"),
  plot.margin = margin(t = -0.3, r = 5.5, b = 5.5, l = 5, unit = "points")
)

# Extract data from the 'total' list
df <- do.call(rbind, lapply(total, function(x) {
  data.frame(
    pre_defined_ci = x$pre_defined_ci,
    power_defined = x$power_defined,
    min_samplesize_for_regressionmodel = x$min_samplesize_for_regressionmodel,
    true_positive_rate = x$true_positive_rate,
    false_positive_rate = x$false_positive_rate,
    true_negative_rate = x$true_negative_rate,
    false_negative_rate = x$false_negative_rate
  )
}))

# Add CI percentiles to df
ci_summary <- sort(unique(df$pre_defined_ci))
df$CI_percentile <- factor(df$pre_defined_ci, levels = ci_summary, labels = c(0.05, 0.25, 0.5, 0.75, 0.95))
file_path <- "Figures/"

# Convert columns to factors for plotting
df$min_samplesize_for_regressionmodel <- as.factor(df$min_samplesize_for_regressionmodel)
df$pre_defined_ci <- as.factor(df$pre_defined_ci)
df$CI_percentile <- as.factor(df$CI_percentile)

# Function to create ROC plot
create_roc_plot <- function(data, ci_percentile, title, show_legend = FALSE) {
  plot <- ggplot(data = data %>% filter(CI_percentile == ci_percentile) %>% arrange(false_positive_rate, true_positive_rate), 
                 aes(x = false_positive_rate, y = true_positive_rate, color = min_samplesize_for_regressionmodel, linetype = min_samplesize_for_regressionmodel), lwd = 2) +
    geom_line(lwd = 1.5) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    xlab("False positive rate") +
    ylab("True positive rate") +
    ylim(0, 1) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(color = "Nmin", linetype = "Nmin") +
    ggtitle(title) +
    scale_color_manual(values = c("10" = colors[1], "20" = colors[2], "30" = colors[3])) +
    scale_linetype_manual(name = "Nmin", values = c("10" = "solid", "20" = "dashed", "30" = "dotted", "BOS" = "solid", "Fixed sample size" = "dashed"))
  
  if (show_legend) {
    plot <- plot + theme(legend.position = "bottom")
  } else {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}

# Create individual ROC plots
g1 <- create_roc_plot(df, 0.05, "5% CIL")
g2 <- create_roc_plot(df, 0.25, "25% CIL")
g3 <- create_roc_plot(df, 0.5, "50% CIL")
g4 <- create_roc_plot(df, 0.75, "75% CIL")
g5 <- create_roc_plot(df, 0.95, "95% CIL", show_legend = TRUE)

# Adjust plot margins
adjust_plot_margins <- function(plot) {
  plot + theme(plot.margin = margin(10, 10, 20, 10))
}

g1 <- adjust_plot_margins(g1)
g2 <- adjust_plot_margins(g2)
g3 <- adjust_plot_margins(g3)
g4 <- adjust_plot_margins(g4) + theme(plot.margin = margin(20, 10, 10, 10))
g5 <- adjust_plot_margins(g5) + theme(plot.margin = margin(20, 10, 10, 10))

# Combine plots into a single plot
combined_plot <- ggarrange(
  g1, g2, g3, g4, g5, 
  labels = c("A", "B", "C", "D", "E"),
  nrow = 2, ncol = 3, 
  common.legend = TRUE, legend = "bottom",
  label.x = 0, 
  label.y = c(0.95, 0.95, 0.95, 0.9, 0.9, 0.9)
)

# Save the combined plot
ggsave(file.path(file_path, "ROCggplot.jpeg"), combined_plot, width = 10, height = 5, dpi = 300)