groups_of_data = 1e5
max_sample_resource = 50
prior_type <- "central_informative_prior" # choose from c("central_informative_prior",
                                          # "central_weakly_informative_prior",
                                          # "offset_weakly_informative_prior",
                                          # "offset_informative_prior",
                                          # "flat_prior")
if (prior_type == "central_informative_prior"){
  prior <- data.frame(
    mu0 <- 0,
    n0 <- 10,
    phi0 <- 1,
    v0 <- 10
  )
}else if (prior_type == "central_weakly_informative_prior"){
  prior <- data.frame(
    mu0 <- 0,
    n0 <- 5,
    phi0 <- 10,
    v0 <- 1
  )
}else if (prior_type == "offset_weakly_informative_prior"){
  prior <- data.frame(
    mu0 <- 1,
    n0 <- 5,
    phi0 <- 10,
    v0 <- 1
  )
}else if (prior_type == "offset_informative_prior"){
  prior <- data.frame(
    mu0 <- 5,
    n0 <- 10,
    phi0 <- 1,
    v0 <- 10
  )
}else if (prior_type == "flat_prior"){
  prior <- data.frame(
    mu0 <- 0,
    n0 <- 1,
    phi0 <- 20,
    v0 <- 1
  )
}
# precreate matrix to store credible interval length
ci_data_all <- matrix(nrow = groups_of_data,ncol = max_sample_resource)
# precreate list to save simulated data
simulated_data_list <- list()
for  (l in 1:groups_of_data){
  set.seed(l)
  simulated_data <- rnorm(max_sample_resource,mean = 0,1)
  simulated_data_list[[l]] <- simulated_data
  for (i in 1:(max_sample_resource)){
    n <- i
    selecteddata <- simulated_data_list[[l]][1:i]
    x_ba <- mean(selecteddata)
    s <- sqrt(1/(n-1) * sum((selecteddata-x_ba)^2))
    # Posterior distribution based on data
    mu1 <- (n0*mu0+n*x_ba)/(n0+n)
    n1 <- n0+n
    v1 <- v0+n
    phi1 <- 1/v1 * ((n-1)*s^2 + v0*phi0 + n*n0/n1 * (x_ba-mu0)^2)
    # 95% credible interval length of the data posterior
    alpha <- 0.05
    df <- n-1
    t_crit <- qt(1 - alpha/2, df)
    ci_data_all[l,i] <- 2*t_crit*sqrt(phi1/n1)
  }
}
pre_defined_ci_percentile <- c(0.05, 0.25, 0.5, 0.75, 0.95)
pre_defined_ci <- sapply(pre_defined_ci_percentile, function(p) quantile(ci_data_all[, max_sample_resource], p))
save(prior,simulated_data_list,ci_data_all,pre_defined_ci,max_sample_resource,pre_defined_ci_percentile, prior_type,
     file = paste("Data/",prior_type,"_statisticstandard",
                  ".RData",sep = ""))
