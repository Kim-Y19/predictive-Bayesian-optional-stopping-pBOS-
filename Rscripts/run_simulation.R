# Load necessary libraries
suppressPackageStartupMessages({
  library(readxl)
  library(gam)
  library(progress)
})

# Source external R script
source("Rscripts/pBOS.R")

# Define file name and load data
file_name <- "central_informative_prior_statisticstandard"
load(file.path("Data", paste0(file_name, ".RData")))

# Define parameters
groups_of_data <- 40
input <- expand.grid(
  pre_defined_ci = as.numeric(pre_defined_ci),
  tolerance_level = seq(0, 1, 0.1),
  Nmin_for_regressionmodel = seq(10, 30, 10)
)

# Initialize result list and progress bar
result <- vector("list", length(input$pre_defined_ci))
pb <- progress_bar$new(
  format = "[:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = length(input$pre_defined_ci),
  clear = FALSE
)

# Main loop to process data
for (i in seq_along(input$pre_defined_ci)) {
  gc()  # Garbage collection to free up memory
  pb$tick()  # Update progress bar
  
  # Call pBOS function with parameters
  res <- pBOS(
    simulated_data_list,  # List of simulated data
    prior,  # Prior information
    Samplesizemax = max_sample_resource,  # Maximum sample size
    pre_defined_ci = input$pre_defined_ci[i],  # Pre-defined credible interval
    groups_of_data = groups_of_data,  # Number of data groups
    MaxDatainfuture = max_sample_resource,  # Maximum data in future
    tolerance_level = input$tolerance_level[i],  # Tolerance level
    Nmin = input$Nmin_for_regressionmodel[i],  # Minimum sample size for regression model
    nSim = 300  # Number of simulations
  )
  
  print(paste("Loop", i, "is finished"))  # Print loop completion message
  result[[i]] <- res  # Save result
}

# Save results for plotting later
save(result, file = file.path("Results", paste0(prior_type, ".RData")))
