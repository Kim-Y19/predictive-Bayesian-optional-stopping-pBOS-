# Load necessary libraries
library(ggplot2)
library(GGally)
library(RColorBrewer)
library(purrr)
library(cowplot)
library(ggpubr)

# Load the data
file_load_name <- "result_20241026_statistics_group400_Regressionyes_prior_var1_offsetmean_0_n0_10_v0_10_centermuTRUE_centervarTRUE_informativemuTRUE_informativevarTRUE_weaklyFALSE_statisticstandard_NEW1025__group_n_1"
load(paste0("Results/", file_load_name, ".RData"))

# Set color palette and line types
colors <- brewer.pal(12, "Paired")
cb_palette <- c("#1B9E77", "#E6AB02", "#666666", "#7570B3", "#D95F02")
colors <- brewer.pal(5, "Dark2")

# Set theme for plots
ptsize <- 16
theme_set(theme_bw())
theme_update(
  axis.text = element_text(size = ptsize, colour = "black", family = "serif"),
  axis.line = element_line(colour = "black", size = 0.25),
  axis.ticks = element_line(colour = "black", size = 0.25),
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
  plot.title = element_text(size = ptsize, colour = "black", family = "serif", face = "plain", hjust = 0.5),
  text = element_text(size = ptsize, colour = "black", family = "serif"),
  panel.background = element_rect(fill = "white", colour = NA),
  plot.background = element_rect(fill = "white", colour = NA),
  legend.background = element_rect(fill = "white", colour = NA)
)

# Initialize variables
total_length <- length(total[[1]]$final_sample_size)
column_names <- c("pre_defined_ci", "power_defined", "min_samplesize_for_regressionmodel", "total_cost_benefit", "ci_percentile", "true_rate", "true_rate2", "propostion", "base_cost_benefit", "run_all_experiments_base", "propostion2", "total_cost_benefit_per_experiment")
benefit_value <- 50  # The benefit for reaching the goal
p = 0 # penalty

# Loop through benefit values
for (cost_index in 1:length(benefit_value)) {
  b <- benefit_value[cost_index]
  cost_benefit_df <- data.frame(matrix(ncol = length(column_names), nrow = length(total)))
  colnames(cost_benefit_df) <- column_names
  pre_defined_ci_all <- sapply(total, function(x) unique(x$pre_defined_ci))
  
  # Loop through total list
  for (loop_index in 1:length(total)) {
    tt <- total[[loop_index]]
    Samplesizemax <- 50
    true_positive <- false_positive <- false_negative <- true_negative <- 0
    true_label <- cost <- benefit <- penalty <- rep(0, length(tt$final_sample_size))
    benefit_lowercost <- rep(0, length(tt$final_sample_size))
    true_sum_true_part <- true_sum_false_part <- 0
    
    # Loop through final sample sizes
    for (i in 1:length(tt$final_sample_size)) {
      if (is.infinite(tt$final_sample_size[i])) {
        if (tt$ci_data_all[i, Samplesizemax] > tt$pre_defined_ci) {
          true_negative <- true_negative + 1
          true_label[i] <- 1
          cost[i] <- tt$stop_sample_size[i]
          benefit_lowercost[i] <- Samplesizemax - tt$stop_sample_size[i]
          penalty[i] <- 0
          if (i <= (total_length / 2)) {
            true_sum_false_part <- true_sum_false_part + 1
          } else {
            true_sum_true_part <- true_sum_true_part + 1
          }
        } else {
          false_negative <- false_negative + 1
          true_label[i] <- 0
          cost[i] <- tt$stop_sample_size[i]
          benefit[i] <- 0
          penalty[i] <- p
        }
      } else {
        if (tt$ci_data_all[i, tt$stop_sample_size[i]] > tt$pre_defined_ci) {
          false_positive <- false_positive + 1
          true_label[i] <- 0
          cost[i] <- Samplesizemax
          benefit[i] <- 0
          penalty[i] <- 0
        } else {
          true_positive <- true_positive + 1
          true_label[i] <- 1
          cost[i] <- tt$stop_sample_size[i]
          benefit[i] <- b
          penalty[i] <- 0
          if (i <= (total_length / 2)) {
            true_sum_false_part <- true_sum_false_part + 1
          } else {
            true_sum_true_part <- true_sum_true_part + 1
          }
        }
      }
    }
    
    # Calculate CI percentiles
    pre_defined_ci <- sort(unique(pre_defined_ci_all))
    if (tt$pre_defined_ci==pre_defined_ci[order(unique(pre_defined_ci))][1]){
      ci_percentile = 0.05
    }else if(tt$pre_defined_ci == pre_defined_ci[order(unique(pre_defined_ci))][2]){
      ci_percentile = 0.25
    }else if(tt$pre_defined_ci == pre_defined_ci[order(unique(pre_defined_ci))][3]){
      ci_percentile = 0.5
    }else if (tt$pre_defined_ci == pre_defined_ci[order(unique(pre_defined_ci))][4]){
      ci_percentile = 0.75
    }else if (tt$pre_defined_ci == pre_defined_ci[order(unique(pre_defined_ci))][5]){
      ci_percentile = 0.95
    }
    
    true_part_cost_benefit <- sum(benefit[1:(total_length / 2)])
    false_part_cost_benefit <- sum(benefit[(total_length / 2 + 1):total_length])
    true_part_cost_trail <- sum(cost[1:(total_length / 2)])
    false_part_cost_trail <- sum(cost[(total_length / 2 + 1):total_length])
    
    total_cost_benefit <- (true_part_cost_benefit * ci_percentile + false_part_cost_benefit * (1 - ci_percentile)) / (false_part_cost_trail * (1 - ci_percentile) + true_part_cost_trail * ci_percentile)
    total_cost_benefit_per_experiment <- ((true_part_cost_benefit + sum(benefit_lowercost[benefit_lowercost])) * (1 - ci_percentile) + (false_part_cost_benefit + sum(benefit_lowercost[1:(total_length / 2)])) * ci_percentile) / (total_length / 2 * ci_percentile + total_length / 2 * (1 - ci_percentile))
    
    cost_benefit_df$total_cost_benefit[loop_index] <- total_cost_benefit
    cost_benefit_df$pre_defined_ci[loop_index] <- tt$pre_defined_ci
    cost_benefit_df$power_defined[loop_index] <- tt$power_defined
    cost_benefit_df$min_samplesize_for_regressionmodel[loop_index] <- tt$min_samplesize_for_regressionmodel
    cost_benefit_df$ci_percentile[loop_index] <- ci_percentile
    cost_benefit_df$total_cost_benefit_per_experiment[loop_index] <- total_cost_benefit_per_experiment
    
    sum_n <- true_negative + false_negative + false_positive + true_positive
    true_positive_rate <- true_positive / (true_positive + false_negative)
    false_postive_rate <- false_positive / (false_positive + true_negative)
    
    cost_benefit_df$true_rate[loop_index] <- (true_sum_true_part * (1 - ci_percentile) + true_sum_false_part * ci_percentile) / total_length * 2
    cost_benefit_df$true_rate2[loop_index] <- (sum(true_label[(1 + total_length / 2):total_length]) * (1 - ci_percentile) + sum(true_label[1:(total_length / 2)]) * ci_percentile) / total_length * 2
    cost_benefit_df$run_all_experiments_base[loop_index] <- (0 * (1 - ci_percentile) / 2 + total_length / 2 * ci_percentile * b) / (Samplesizemax * total_length / 2 * ci_percentile + Samplesizemax * total_length / 2 * (1 - ci_percentile))
  }
  
  # Calculate baseline cost-benefit
  gg_count <- 0
  get_baseline_bayesianstop <- numeric(length(unique(cost_benefit_df$pre_defined_ci)))
  for (i in unique(cost_benefit_df$pre_defined_ci)) {
    gg_count <- gg_count + 1
    get_baseline_bayesianstop[gg_count] <- cost_benefit_df[cost_benefit_df$pre_defined_ci == i & cost_benefit_df$power_defined == 0.0,]$total_cost_benefit[1]
    cost_benefit_df[cost_benefit_df$pre_defined_ci == i,]$base_cost_benefit <- get_baseline_bayesianstop[gg_count]
  }
  
  # Calculate proportions
  for (i in 1:length(cost_benefit_df$pre_defined_ci)) {
    cost_benefit_df$propostion[i] <- (cost_benefit_df$total_cost_benefit[i] - cost_benefit_df$base_cost_benefit[i]) / cost_benefit_df$base_cost_benefit[i]
    cost_benefit_df$propostion2[i] <- (cost_benefit_df$run_all_experiments_base[i] - cost_benefit_df$base_cost_benefit[i]) / cost_benefit_df$base_cost_benefit[i]
  }
  
  # Plotting
  group_Nmin <- unique(cost_benefit_df$min_samplesize_for_regressionmodel)
  cost_benefit_df$min_samplesize_for_regressionmodel <- factor(cost_benefit_df$min_samplesize_for_regressionmodel, levels = group_Nmin)
  
  plot_cost_benefit <- function(data, ci_percentile, title) {
    ggplot(data = data[data$ci_percentile == ci_percentile,], aes(x = power_defined, y = propostion, color = min_samplesize_for_regressionmodel, linetype = min_samplesize_for_regressionmodel)) +
      geom_line(aes(linetype = min_samplesize_for_regressionmodel, color = min_samplesize_for_regressionmodel), lwd = 1.5) +
      geom_hline(aes(yintercept = 0, linetype = "BOS", color = "BOS"), lwd = 1.5, show.legend = TRUE) +
      geom_hline(aes(yintercept = unique(data[data$ci_percentile == ci_percentile,]$propostion2), linetype = "Fixed sample size", color = "Fixed sample size"), lwd = 1.5, show.legend = TRUE) +
      scale_color_manual(name = "Nmin", values = c("10" = colors[1], "20" = colors[2], "30" = colors[3], "BOS" = "gray", "Fixed sample size" = "gray")) +
      scale_linetype_manual(name = "Nmin", values = c("10" = "solid", "20" = "dashed", "30" = "dotted", "BOS" = "solid", "Fixed sample size" = "dashed")) +
      xlab("Tolerance level") +
      ylab("Proportion increase") +
      ylim(-1, 1) +
      ggtitle(title) +
      labs(color = "Nmin", linetype = "Nmin")
  }
  
  g1 <- plot_cost_benefit(cost_benefit_df, 0.05, "5% CIL")
  g2 <- plot_cost_benefit(cost_benefit_df, 0.25, "25% CIL")
  g3 <- plot_cost_benefit(cost_benefit_df, 0.5, "50% CIL")
  g4 <- plot_cost_benefit(cost_benefit_df, 0.75, "75% CIL")
  g5 <- plot_cost_benefit(cost_benefit_df, 0.95, "95% CIL") +
    guides(color = guide_legend(nrow = 3, title.position = "top", override.aes = list(size = 1.1)), linetype = guide_legend(override.aes = list(size = 1.1))) +
    theme(legend.position = "bottom")
  
  # Extract the legend
  leg <- get_legend(g5)
  
  # Convert to a ggplot and print
  g6 <- as_ggplot(leg)
  
  # Combine the plots with transparent background
  combined_plot <- ggarrange(
    g1, g2, g3, g6, g4, g5,
    labels = c("A", "B", "C", "", "D", "E"), font.label = list(size = ptsize, color = "white"),
    nrow = 2, ncol = 3,
    common.legend = TRUE, legend = "none",
    label.x = 0,
    label.y = c(0.95, 0.95, 0.95, 0.9, 0.9, 0.9)
  ) + theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )
  
  # Save the plot with ggsave
  ggsave(paste("Figures/Cost_benefit_per_trial_fig", cost_index, ".png",sep = ""), combined_plot, width = 11, height = 6, dpi = 300)
}