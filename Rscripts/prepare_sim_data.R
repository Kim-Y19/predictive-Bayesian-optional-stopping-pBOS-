# Set parameters
groups_of_data <- 1e6  # Number of groups of data to simulate
max_sample_resource <- 50  # Maximum sample size for each group
prior_type <- "central_informative_prior"  # Type of prior to use
# Define prior based on prior_type
prior <- switch(prior_type,
                "central_informative_prior" = data.frame(mu0 = 0, n0 = 10, phi0 = 1, v0 = 10),
                "central_weakly_informative_prior" = data.frame(mu0 = 0, n0 = 5, phi0 = 10, v0 = 1),
                "offset_weakly_informative_prior" = data.frame(mu0 = 1, n0 = 5, phi0 = 10, v0 = 1),
                "offset_informative_prior" = data.frame(mu0 = 5, n0 = 10, phi0 = 1, v0 = 10),
                "flat_prior" = data.frame(mu0 = 0, n0 = 1, phi0 = 20, v0 = 1)
)

# Pre-create matrix to store credible interval lengths
ci_data_all <- matrix(nrow = groups_of_data, ncol = max_sample_resource)

# Pre-create list to save simulated data
simulated_data_list <- vector("list", groups_of_data)

# Generate simulated data and calculate credible intervals
for (l in 1:groups_of_data) {
  set.seed(l)  # Set seed for reproducibility
  simulated_data <- rnorm(max_sample_resource, mean = 0, sd = 1)  # Generate random data
  simulated_data_list[[l]] <- simulated_data  # Save simulated data
  
  for (i in 1:max_sample_resource) {
    n <- i  # Current sample size
    selecteddata <- simulated_data[1:n]  # Select data up to current sample size
    x_ba <- mean(selecteddata)  # Calculate sample mean
    s <- sqrt(1 / (n - 1) * sum((selecteddata - x_ba)^2))  # Calculate sample standard deviation
    
    # Posterior distribution based on data
    mu1 <- (prior$n0 * prior$mu0 + n * x_ba) / (prior$n0 + n)  # Update mean posterior distribution
    n1 <- prior$n0 + n  # Update sample size in posterior distribution
    v1 <- prior$v0 + n  
    phi1 <- 1 / v1 * ((n - 1) * s^2 + prior$v0 * prior$phi0 + n * prior$n0 / n1 * (x_ba - prior$mu0)^2)  # Update variance
    
    # 95% credible interval length of the data posterior
    alpha <- 0.05  # Significance level
    df <- n - 1  # Degrees of freedom
    t_crit <- qt(1 - alpha / 2, df)  # Critical value from t-distribution
    ci_data_all[l, i] <- 2 * t_crit * sqrt(phi1 / n1)  # Calculate credible interval length
  }
}

# Define percentiles for pre-defined credible intervals
pre_defined_ci_percentile <- c(0.05, 0.25, 0.5, 0.75, 0.95)  # Percentiles to calculate
pre_defined_ci <- sapply(pre_defined_ci_percentile, function(p) quantile(ci_data_all[, max_sample_resource], p))  # Calculate percentiles

# Save results
save(prior, simulated_data_list, ci_data_all, pre_defined_ci, max_sample_resource, pre_defined_ci_percentile, prior_type,
     file = paste0("Data/", prior_type, "_statisticstandard.RData"))  # Save data to file
