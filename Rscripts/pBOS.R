pBOS <- function(
    simulated_data_list, # List of simulated data
    prior, # Prior information
    Samplesizemax = 500, # Maximum sample size
    pre_defined_ci = 1,  # Pre-defined credible interval
    groups_of_data = 100,  # Number of data groups
    MaxDatainfuture = 100,  # Maximum data in future
    tolerance_level = 1,  # Tolerance level
    Nmin = 10,  # Minimum sample size for pBOS
    min_samplesize_for_regressionmodel = 3, #minimum sample size for regression model input
    nSim = 1000  # Number of simulations
) {
  # Extract prior parameters
  mu0 <- prior$mu0
  n0 <- prior$n0
  phi0 <- prior$phi0
  v0 <- prior$v0
  
  # Set seed for reproducibility
  set.seed(1)
  # mean_theta = v0/2*(2/(v0*phi0))
  # theta <- rgamma(10e5, shape = v0/2, scale  = 2/(v0*phi0))  # Generate gamma-distributed values
  # mu <- rnorm(10e5, mean = mu0, sd = sqrt(1/(n0* theta)))  # Generate normal-distributed values
  # phi = 1/theta
  
  # Initialize matrices and lists for storing results
  proportion_ci_pred_reaches <- matrix(nrow = groups_of_data, ncol = Samplesizemax)
  final_sample_size <- c()
  stop_sample_size <- c()
  ci_predited_withpower_at_i_samplesize <- matrix(nrow = groups_of_data, ncol = MaxDatainfuture)
  ci_data <- matrix(nrow = groups_of_data, ncol = MaxDatainfuture)
  false_label <- c()
  true_label <- c()
  
  # Classify data groups based on predefined CI
  for (l in 1:length(simulated_data_list)) {
    if(any(!is.na(ci_data_all[l,]) & ci_data_all[l,] < pre_defined_ci )){
      false_label[l] <- l
    } else {
      true_label[l] <- l
    }
  }
  
  ratio_generated_data_variance_original_data_variance <- matrix(nrow = groups_of_data, ncol = MaxDatainfuture)
  n_select <- groups_of_data / 2

  # Randomly select true and false labels
  true_select <- sample(na.omit(true_label), n_select)
  false_select <- sample(na.omit(false_label), n_select)
  combine_true_false <- c(false_select, true_select)
  sum_selected_data_list <- vector("list", groups_of_data)
  ci_data_all <- matrix(nrow = groups_of_data, ncol = MaxDatainfuture)
  # Main loop to process each data group
  for (l in 1:groups_of_data) {
    set.seed(l)
    sum_selected_data_list[[l]] <- simulated_data_list[[combine_true_false[l]]]
    realdata <- simulated_data_list[[combine_true_false[l]]]
    data_to_generate <- vector("list", nSim)
    theta_g = c()
    mu_g = c()
    phi_g = c()
    xg_ba <- matrix(nrow = MaxDatainfuture, ncol = nSim)
    sg <- matrix(nrow = MaxDatainfuture, ncol = nSim)
    mug_post <- matrix(nrow = MaxDatainfuture, ncol = nSim)
    phig_post <- matrix(nrow = MaxDatainfuture, ncol = nSim)
    ci_g_post <- matrix(nrow = MaxDatainfuture, ncol = nSim)
    
    # Loop to generate data and calculate CI
    for (i in 1:Samplesizemax) {
      n <- i
      selecteddata <- realdata[1:n]
      x_ba <- mean(selecteddata)
      s <- sqrt(1 / (n - 1) * sum((selecteddata - x_ba)^2))
      
      # Calculate posterior parameters
      mu1 <- (n0 * mu0 + n * x_ba) / (n0 + n)
      n1 <- n0 + n
      v1 <- v0 + n
      phi1 <- 1 / v1 * ((n - 1) * s^2 + v0 * phi0 + n * n0 / n1 * (x_ba - mu0)^2)
      
      # Calculate CI length
      alpha <- 0.05
      df <- n - 1
      t_crit <- qt(1 - alpha / 2, df)
      ci_data[l, i] <- 2 * t_crit * sqrt(phi1 / n1)
      ratio_generated_data_variance_original_data_variance[l, i] <- phi1 / s^2
      
      # Check if CI length meets the predefined threshold
      if (!is.na(ci_data[l, i]) & ci_data[l, i] <= pre_defined_ci) {
        print("data has reached the precision target")
        final_sample_size[l] <- i
        stop_sample_size[l] <- i
        break
      } else {
        theta_g <- rgamma(nSim, shape = v1 / 2, scale = 2 / (v1 * phi1))
        for (k in 1:nSim) {
          set.seed(k)
          mu_g[k] <- rnorm(1, mean = mu1, sd = sqrt(1 / (n1 * theta_g[k])))
          phi_g[k] <- 1 / theta_g[k]
          for (j in 1:MaxDatainfuture) {
            data_to_generate[[k]] <- rnorm(j, mean = mu_g[k], sd = sqrt(phi_g[k]))
            xg_ba[j, k] <- mean(data_to_generate[[k]])
            sg[j, k] <- sqrt(1 / (j - 1) * sum((data_to_generate[[k]] - xg_ba[j, k])^2))
            mug_post[j, k] <- (n0 * mu0 + j * xg_ba[j, k]) / (n0 + j)
            ng_post <- n0 + j
            vg_post <- v0 + j
            phig_post[j, k] <- 1 / vg_post * ((j - 1) * sg[j, k]^2 + v0 * phi0 + j * n0 / ng_post * (xg_ba[j, k] - mu0)^2)
            ci_g_post[j, k] <- 2 * t_crit * sqrt(phig_post[j, k] / ng_post)
          }
        }
        
        # Regression model
        if (i >= Nmin) {
          percentile_offset <- quantile(ci_g_post[MaxDatainfuture, ], tolerance_level) - median(ci_g_post[MaxDatainfuture, ])
          real_data_ci <- c()
          num_real_data <- c()
          num_predicted_data <- c()
          median_predicted_ci <- c()
          percentile_predicted_ci <- c()
          g <- 1
          for (p in 1:i) {
            for (j in 1:MaxDatainfuture) {
              real_data_ci[g] <- ifelse(j <= i, ci_data[l, j], NA)
              if (anyNA(ci_g_post[j, ])) {
                percentile_predicted_ci[g] <- NA
                median_predicted_ci[g] <- NA
              } else {
                percentile_predicted_ci[g] <- quantile(ci_g_post[j, ], tolerance_level, na.rm = TRUE)
                median_predicted_ci[g] <- median(ci_g_post[j, ])
              }
              num_real_data[g] <- p
              num_predicted_data[g] <- j
              g <- g + 1
            }
          }
          
          all_data <- data.frame(
            y = 1 / real_data_ci^2,
            x1 = 1 / median_predicted_ci^2,
            x2 = num_real_data,
            x3 = num_predicted_data
          )
          data <- all_data[complete.cases(all_data), ]
          sorted_data <- data[order(data$x2, data$x3), ]
          train_group <- sorted_data[(sorted_data$x2 <= i & 
                                        sorted_data$x3 >= min_samplesize_for_regressionmodel & 
                                        sorted_data$x3 <= i & 
                                        sorted_data$x2 >= min_samplesize_for_regressionmodel), ]
          test_group <- all_data[all_data$x2 == i, ]
          
          # Initialize an empty list to store the predictions
          model <- lm(y ~ x1 + x2 + x3, data = train_group)
          predictions = predict(model,newdata=test_group)
          test_group$predictions = predictions
          # Adjust CI based on predictions
          # if (length(test_group[test_group$x3 == MaxDatainfuture, ]$predictions) > 0 && 
          #     !is.na(test_group[test_group$x3 == MaxDatainfuture, ]$predictions) && 
          #     test_group[test_group$x3 == MaxDatainfuture, ]$predictions > 0) 
          if(test_group[test_group$x3 == MaxDatainfuture,]$predictions>0)
            {
            ecdf_function <- ecdf(ci_g_post[MaxDatainfuture, ] + 
                                    sqrt(1 / test_group[test_group$x3 == MaxDatainfuture, ]$predictions) - 
                                    median(ci_g_post[MaxDatainfuture, ]))
            ci_predited_withpower_at_i_samplesize[l, i] <- sqrt(1 / test_group[test_group$x3 == MaxDatainfuture, ]$predictions) + 
              percentile_offset
          } else {
            ecdf_function <- ecdf(ci_g_post[MaxDatainfuture, ])
            ci_predited_withpower_at_i_samplesize[l, i] <- quantile(ci_g_post[MaxDatainfuture, ], tolerance_level)
          }
          # Adjust CI based on tolerance level
          if (tolerance_level == 0) {
            ci_predited_withpower_at_i_samplesize[l, i] <- 0
          } else if (tolerance_level == 1) {
            ci_predited_withpower_at_i_samplesize[l, i] <- 1000
          } 
          # else {
          #   ci_predited_withpower_at_i_samplesize[l, i] <- NA
          # }
          
          # Check if experiments have reached the initial decision sample size
          if (i < Nmin) {
            # Experiments have not reached the initial pBOS stop decision sample size
          } else if (!is.na(ci_predited_withpower_at_i_samplesize[l, i]) & ci_predited_withpower_at_i_samplesize[l, i] > pre_defined_ci) {
            final_sample_size[l] <- Inf
            stop_sample_size[l] <- i
            break
          } else if (i == Samplesizemax) {
            # browser()
            final_sample_size[l] <- i
            stop_sample_size[l] <- i
          }
        }
      }
    }
    print(paste("Now it is the group", l))
  }
  
  # Calculate CI for all data groups
  for (l in 1:groups_of_data) {
    for (i in 1:Samplesizemax) {
      n <- i
      selecteddata <- sum_selected_data_list[[l]][1:n]
      x_ba <- mean(selecteddata)
      s <- sqrt(1 / (n - 1) * sum((selecteddata - x_ba)^2))
      
      mu1 <- (n0 * mu0 + n * x_ba) / (n0 + n)
      n1 <- n0 + n
      v1 <- v0 + n
      phi1 <- 1 / v1 * ((n - 1) * s^2 + v0 * phi0 + n * n0 / n1 * (x_ba - mu0)^2)
      
      alpha <- 0.05
      df <- n - 1
      t_crit <- qt(1 - alpha / 2, df)
      ci_data_all[l, i] <- 2 * t_crit * sqrt(phi1 / n1)
      
    }
  }

  # Calculate true and false positives/negatives
  false_positive <- 0
  true_positive <- 0
  true_negative <- 0
  false_negative <- 0
  true_label <- c()

  for (i in 1:length(final_sample_size)) {
    if (is.infinite(final_sample_size[i])) {
      if (ci_data_all[i, Samplesizemax] > pre_defined_ci) {
        true_negative <- true_negative + 1
        true_label[i] <- 1
      } else {
        false_negative <- false_negative + 1
        true_label[i] <- 0
      }
    } else {
      if (ci_data_all[i, stop_sample_size[i]] > pre_defined_ci) {
        false_positive <- false_positive + 1
        true_label[i] <- 0
      } else {
        true_positive <- true_positive + 1
        true_label[i] <- 1
      }
    }
  }

  false_positive_rate <- false_positive / (false_positive + true_negative)
  false_negative_rate <- false_negative / length(final_sample_size)
  true_positive_rate <- true_positive / (true_positive + false_negative)
  true_negative_rate <- true_negative / length(final_sample_size)
  # browser()
  # Return results
  results <- list(
    Nmin = Nmin,
    pre_defined_ci = pre_defined_ci,
    tolerance_level = tolerance_level,
    groups_of_data = groups_of_data,
    ci_data_all = ci_data_all,
    final_sample_size = final_sample_size,
    stop_sample_size = stop_sample_size,
    true_positive_rate = true_positive_rate,
    false_positive_rate = false_positive_rate,
    false_negative_rate = false_negative_rate,
    true_negative_rate = true_negative_rate
  )
  return(results)
}
