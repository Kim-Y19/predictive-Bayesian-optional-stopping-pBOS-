BOS <- function(data,
                experiment_data_groups,
                prior,
                Regression = TRUE,
                regression_model_choice =1,
                Samplesizemax = 500,
                pre_defined_ci = 1,
                groups_of_data = 100,
                MaxAddDataInFuture = 100,
                power_defined = 1,
                number_of_initial_sample = 10,
                min_samplesize_for_regressionmodel = 100,
                nSim = 1000,
                n_picks = 1e5,
                gg = 0,
                file_path) {
  mu0 <- prior$mu0
  n0 <- prior$n0
  phi0 <- prior$phi0
  v0 <- prior$v0
  all_row_data_math_median_data_ci = list()
  all_row_data_math_percentile_data_ci = list()

  # Need to check here again for gamma distribution
  set.seed(1)
  mean_theta = v0/2*(2/(v0*phi0))
  theta <- rgamma(n_picks, shape = v0/2, scale  = 2/(v0*phi0))  # Generate gamma-distributed values
  mu <- rnorm(n_picks, mean = mu0, sd = sqrt(1/(n0* theta)))  # Generate normal-distributed values
  phi = 1/theta

  proportion_ci_pred_reaches = matrix(nrow = groups_of_data, ncol = Samplesizemax)
  realdata_list = list()
  final_sample_size = c()
  stop_sample_size = c()
  rsquared = matrix(nrow = groups_of_data, ncol = MaxAddDataInFuture)
  check_percentile = matrix(nrow = groups_of_data, ncol = MaxAddDataInFuture)
  ci_predited_withpower_at_i_samplesize = matrix(nrow = groups_of_data, ncol = MaxAddDataInFuture)
  predicted_median_summary = matrix(nrow = groups_of_data, ncol = MaxAddDataInFuture)
  without_regression_median = matrix(nrow = groups_of_data, ncol = MaxAddDataInFuture)
  ci_data = matrix(nrow = groups_of_data, ncol = MaxAddDataInFuture)
  true_group_data = list()
  false_group_data = list()
  false_label = c()
  true_label = c()
  for(l in 1:length(experiment_data_groups)){
    if(any(!is.na(ci_data_all[l,]) & ci_data_all[l,] < pre_defined_ci )){
      false_label[l] = l
    }else{
      true_label[l] = l
    }
  }
  ci_data_all = matrix(nrow = groups_of_data, ncol = MaxAddDataInFuture)
  ratio_generated_data_variance_original_data_variance= matrix(nrow = groups_of_data, ncol = MaxAddDataInFuture)
  n_select =groups_of_data/2
  # Randomly select some elements from the vector
  # browser()
  true_select <- sample(na.omit(true_label), n_select )
  false_select <- sample(na.omit(false_label), n_select )
  combine_true_false = c(false_select,true_select)
  sum_selected_data_list = list()
  for  (l in 1:groups_of_data){
    set.seed(l)
    sum_selected_data_list[[l]] = experiment_data_groups[[combine_true_false[l]]]
    realdata = experiment_data_groups[[combine_true_false[l]]]
    data_to_generate = list()
    theta_g = c()
    mu_g = c()
    phi_g = c()
    xg_ba = matrix(nrow = MaxAddDataInFuture, ncol = nSim)
    sg = matrix(nrow = MaxAddDataInFuture, ncol = nSim)
    
    mug_post = matrix(nrow = MaxAddDataInFuture, ncol = nSim)
    phig_post = matrix(nrow = MaxAddDataInFuture, ncol = nSim)
    ci_g_post = matrix(nrow = MaxAddDataInFuture, ncol = nSim)
    
    all_simulated_mu = list()
    all_simulated_phi = list()
    
    all_post_mu = list()
    all_post_phi = list()
    all_sg = list()
    all_ci_g_post = list()
    all_generated_in_ci = list()
    # Plot the generated values
    for (i in 1:(Samplesizemax)){
      n = i
      selecteddata = realdata[1:n]
      realdata_list[[i]] = selecteddata
      x_ba = mean(selecteddata)
      s = sqrt(1/(n-1) * sum((selecteddata-x_ba)^2))
      # Data posterior is
      mu1 = (n0*mu0+n*x_ba)/(n0+n)
      n1 = n0+n
      v1 = v0+n
      phi1 = 1/v1 * ((n-1)*s^2 + v0*phi0 + n*n0/n1 * (x_ba-mu0)^2)
      # ci_length of the data posterior is 
      alpha <- 0.05
      df <- n-1
      t_crit <- qt(1 - alpha/2, df)
      ci_data[l,i] = 2*t_crit*sqrt(phi1/n1)
      # ratio_generated_data_variance_original_data_variance[l,i] = 
      #   (n-1)*s/(n + v0) +
      #   (v0*phi0 + n*n0/(n +n0)*(xbar - phi0))/s
      ratio_generated_data_variance_original_data_variance[l,i] = phi1/s^2
      # Here is the first check, wehther the collected data could provide enough precision
      if (!is.na(ci_data[l,i]) & ci_data[l,i]<= pre_defined_ci){
        print("data has reached the precision target")
        final_sample_size[l] = i
        stop_sample_size[l] =  i
        break
      }
      else{
        theta_g <- rgamma(nSim, shape = v1/2, scale = 2/(v1*phi1))  # Generate gamma-distributed values
        for(k in 1:nSim){ # k is how many groups of simulated samples in future
          set.seed(k)
          mu_g[k] <- rnorm(1, mean = mu1, sd = sqrt(1/(n1* theta_g[k])))  # Generate normal-distributed values
          phi_g[k] = 1/theta_g[k]
          for (j in 1:MaxAddDataInFuture){
            data_to_generate[[k]] = rnorm(j, mean = mu_g[k], sd = sqrt(phi_g[k]))
            xg_ba[j,k] =  mean(data_to_generate[[k]])
            sg[j,k] = sqrt(1/(j-1) * sum((data_to_generate[[k]]-xg_ba[j,k])^2))
            mug_post[j,k] = (n0*mu0+j*xg_ba[j,k])/(n0+j)
            ng_post = n0+j
            vg_post = v0+j
            phig_post[j,k] = 1/vg_post * ((j-1)*sg[j,k]^2 + v0*phi0 + j*n0/ng_post * (xg_ba[j,k]-mu0)^2)
            ci_g_post[j,k] = 2*t_crit*sqrt(phig_post[j,k]/ng_post)
          }
        }
        # Regression model
        if (i>=min_samplesize_for_regressionmodel){
          if(Regression == TRUE){
            # Offset at max sample size between the ci predicted with a specific power and the median of the predicted ci
            # browser()
            percentile_offset = quantile(ci_g_post[MaxAddDataInFuture,],power_defined) - median(ci_g_post[MaxAddDataInFuture,])
            real_data_ci = c()
            predicted_ci = c()
            num_real_data = c()
            num_predicted_data = c()
            median_predicted_ci = c()
            percentile_predicted_ci = c()
            # @ Xiaomi check once more here
            g = 1
              for (p in 1:i){
                for (j in 1:MaxAddDataInFuture){
                  if(j<=i){
                    real_data_ci[g] = ci_data[l,j]
                  }else {
                    real_data_ci[g] = NA
                  }
                  if(anyNA(ci_g_post[j,])){
                    percentile_predicted_ci[g] = NA
                    median_predicted_ci[g] = NA
                  }else{
                    # browser()
                    percentile_predicted_ci[g] = quantile(ci_g_post[j,],power_defined,na.rm = TRUE)
                    median_predicted_ci[g] = median(ci_g_post[j,])
                  }
                  num_real_data[g] = p
                  num_predicted_data[g] = j
                  g = g+1
                }
              }
            
              all_data = data.frame(y = 1/real_data_ci^2,y2 = 1/real_data_ci,x1 = 1/percentile_predicted_ci^2, x2 = num_real_data,
                                    x3 = num_predicted_data, x4 = 1/median_predicted_ci^2,x5 = power_defined,
                                    x6 = percentile_predicted_ci, x7 = median_predicted_ci,x8 = 1/median_predicted_ci)
              data <- all_data[complete.cases(all_data), ]
              sorted_data <- data[order(data$x2, data$x3), ]
              train_group = sorted_data[(sorted_data$x2<=i & 
                                           sorted_data$x3 >= number_of_initial_sample & 
                                           sorted_data$x3 <= i & 
                                           sorted_data$x2>= number_of_initial_sample),]
              test_group = all_data[all_data$x2 == i,]
              ################
              # Initialize an empty list to store the predictions
              # Loop through each iteration
              # Generate new data for each iteration
              if(regression_model_choice ==1)
              {model <- lm(y ~ x2 + x3 + x4, data = train_group)
              }else{
                # Fit the model using only significant predictors
                model <- lm(y ~ 1, data = train_group)  # Start with an intercept-only model
                while (TRUE) {
                  previous_model <- model  # Store the previous model
                  # Perform forward selection and add the most significant predictor
                  model <- step(model, scope = formula(~ .  + x2 + x3 + x4), direction = "forward")
                  # Exit the loop if no additional predictor was added
                  if (all(names(model$coefficients) %in% names(previous_model$coefficients))) {
                    break
                  }
                }
              }
              
              rsquared[l,i] = summary(model)$r.squared
              # Create the empirical cumulative distribution function
              # browser()
              predictions = predict(model,newdata=test_group)
              test_group$predictions = predictions
              # print(test_group[test_group$x3 == MaxAddDataInFuture,])
              # print(sqrt(1/test_group[test_group$x3 == MaxAddDataInFuture,]$predictions))
              if(test_group[test_group$x3 == MaxAddDataInFuture,]$predictions>0){
                ecdf_function <- ecdf(ci_g_post[MaxAddDataInFuture,]+
                                        sqrt(1/test_group[test_group$x3 == MaxAddDataInFuture,]$predictions) -
                                        median(ci_g_post[MaxAddDataInFuture,]))
                
                ci_predited_withpower_at_i_samplesize[l,i] = sqrt(1/test_group[test_group$x3 == MaxAddDataInFuture,]$predictions) + 
                  percentile_offset
              }else {
                ecdf_function <- ecdf(ci_g_post[MaxAddDataInFuture,])
                ci_predited_withpower_at_i_samplesize[l,i] = quantile(ci_g_post[MaxAddDataInFuture,],power_defined)
              }
              check_percentile[l,i] = ecdf_function(pre_defined_ci)
              # browser()
              
              
              # if(power_defined ==0){
              #   ci_predited_withpower_at_i_samplesize[l,i] = 0
              # }else if(power_defined == 1){
              #   ci_predited_withpower_at_i_samplesize[l,i] = 1000
              # }
              # This two varaible is a bit weird. @ Xiaomi, not sure what exactly it is used later.
              # predicted_median_summary[l,i] = sqrt(1/test_group[test_group$x3 == MaxAddDataInFuture,]$predictions)
              # without_regression_median[l,i] = median(ci_g_post[MaxAddDataInFuture,])
            }
          else {
            browser()
            ci_predited_withpower_at_i_samplesize[l,i] = quantile(ci_g_post[MaxAddDataInFuture,],power_defined)}
          if(power_defined ==0){
              ci_predited_withpower_at_i_samplesize[l,i] = 0
          }else if(power_defined == 1){
              ci_predited_withpower_at_i_samplesize[l,i] = 1000
          }
        
        }else {
            ci_predited_withpower_at_i_samplesize[l,i] = NA
        }
        # Have a look at the logic again!@Xiaomi
        if (i<min_samplesize_for_regressionmodel){
          # print("experiments have not reached the initial decision sample size")
        } else if(!is.na(ci_predited_withpower_at_i_samplesize[l,i]) &  ci_predited_withpower_at_i_samplesize[l,i] > pre_defined_ci){
          # print("experiments should end early as the goal is not likely to be reached by the maximum sample to be collected")
          final_sample_size[l] =  Inf
          stop_sample_size[l] =  i
          break
        }else if(i == Samplesizemax){
          # Maybe set a broswer() here
          final_sample_size[l] =  i
          stop_sample_size[l] =  i
        }
        all_simulated_mu[[i]] = mu_g
        all_simulated_phi[[i]] = phi_g
        all_ci_g_post[[i]] = ci_g_post
        all_post_mu[[i]] = mug_post
        all_post_phi[[i]] = phig_post
        all_sg[[i]] = sg
      }
    }
    print(paste("now it is the group ",l))
    # browser()
  }

  for  (l in 1:groups_of_data){
    for (i in 1:(Samplesizemax)){
      n = i
      selecteddata = sum_selected_data_list[[l]][1:n]
      x_ba = mean(selecteddata)
      s = sqrt(1/(n-1) * sum((selecteddata-x_ba)^2))
      # Data posterior is
      mu1 = (n0*mu0+n*x_ba)/(n0+n)
      n1 = n0+n
      v1 = v0+n
      phi1 = 1/v1 * ((n-1)*s^2 + v0*phi0 + n*n0/n1 * (x_ba-mu0)^2)
      # ci_length of the data posterior is 
      alpha <- 0.05
      
      df <- n-1
      t_crit <- qt(1 - alpha/2, df)
      ci_data_all[l,i] = 2*t_crit*sqrt(phi1/n1)
    }
  }
  # stop the experiment early is negative
  # not stop the experiment early is positive
  false_positive = 0
  true_positive = 0
  true_negative = 0 
  false_negative = 0
  true_label = c()
  for (i in 1:length(final_sample_size)){
    if (is.infinite(final_sample_size[i])){
      if (ci_data_all[i,Samplesizemax] > pre_defined_ci){
        true_negative = true_negative + 1
        true_label[i] = 1
      }else{
        false_negative = false_negative + 1
        true_label[i] = 0
      }
    } else{
      if (ci_data_all[i,stop_sample_size[i]] > pre_defined_ci)
      {
        false_positive = false_positive + 1
        true_label[i] = 0
      }else{
        true_positive = true_positive + 1
        true_label[i] = 1
      }
    }
  }
  false_positive_rate = false_positive/(false_positive + true_negative)
  false_negative_rate = false_negative/length(final_sample_size)
  true_positive_rate = true_positive/(true_positive + false_negative)
  
  true_negative_rate = true_negative/length(final_sample_size)
  true_rate = (true_positive + true_negative)/length(final_sample_size)
  print(paste("True rate is:",true_rate))
  
# plots  
  jpeg(filename = file_path, width = 10, height = 10,units = "cm",res = 300) 
  plot(NULL, xlim = c(0, 50), ylim = c(0,2), xlab = 'Sample size', ylab = 'CI length', 
       main = paste("ci_thre =",pre_defined_ci,
                    ", power = ", power_defined,
                    "\ninitial stop sample = ", min_samplesize_for_regressionmodel,
                    ", correct rate = ",true_rate),
       cex.main = 0.9)
  # Loop through each row of the matrix
  for (i in 1:nrow(ci_data_all)) {
    x_values <- seq(0, 50, length.out = ncol(ci_data_all))
    y_values <- ci_data_all[i, ]
    
    # Split x and y values before and after the threshold value t[i]
    x_values_solid <- x_values[x_values <= stop_sample_size[i]]
    y_values_solid <- y_values[x_values <= stop_sample_size[i]]
    x_values_dashed <- x_values[x_values > stop_sample_size[i]]
    y_values_dashed <- y_values[x_values > stop_sample_size[i]]
    
    # Plot the line with the appropriate style
    lines(x_values_solid, y_values_solid, type = 'l', col = i, lty = 'solid',lwd = 3)
    lines(x_values_dashed, y_values_dashed, type = 'l', col = i, lty = 'dashed')
    # Add points at the threshold t[i]
    
  }
  # Add a horizontal line at y = 0.3
  abline(h = pre_defined_ci, col = 'red')
  dev.off()  # Close the graphics device to save the file
  results = list(
                 final_sample_size = final_sample_size,
                 proportion_ci_pred_reaches = proportion_ci_pred_reaches,
                 pre_defined_ci = pre_defined_ci,
                 ci_data_all = ci_data_all,
                 stop_sample_size = stop_sample_size,
                 true_rate = true_rate,
                 groups_of_data = groups_of_data,
                 min_samplesize_for_regressionmodel = min_samplesize_for_regressionmodel,
                 power_defined = power_defined,
                 true_positive_rate = true_positive_rate,
                 false_positive_rate = false_positive_rate,
                 false_negative_rate = false_negative_rate,
                 true_negative_rate = true_negative_rate,
                ci_predited_withpower_at_i_samplesize = ci_predited_withpower_at_i_samplesize,
                rsquared = rsquared,
                check_percentile = check_percentile,
                only_truepositive = true_positive/length(final_sample_size),
                only_falsepositive = false_positive/length(final_sample_size),
                true_label = true_label,
                ratio_generated_data_variance_original_data_variance = ratio_generated_data_variance_original_data_variance)
}
