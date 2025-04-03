# Load necessary libraries
library("readxl")
library("gam")
library("progress")

# Source external R script
source("Rscripts/pBOS.R")

# Define file name and load data
file_name <- "central_informative_prior_statisticstandard"
load(paste("Data/", file_name, ".RData", sep = ""))

# Define parameters
groups_of_data <- 400
input <- expand.grid(
  pre_defined_ci = as.numeric(pre_defined_ci),
  tolerance_level = seq(0, 1, 0.1),
  Nmin_for_regressionmodel = seq(10, 30, 10)
)

# Initialize result list and progress bar
result <- replicate(length(input$pre_defined_ci), data.frame())
index <- 0
pb <- progress_bar$new(
  format = "[:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = length(input$pre_defined_ci),
  clear = FALSE
)

# Main loop to process data
for (i in 1:length(input$pre_defined_ci)) {
  gc()  # Garbage collection to free up memory
  index <- index + 1
  pb$tick()  # Update progress bar
  
  # Define file path for saving figures
  # file_path <- file.path(paste("Figures/", file_name, sep = ""), paste0("Figure", index, ".jpeg"))
  
  # Call pBOS function with parameters
  res <- pBOS(
    simulated_data_list,  # List of simulated data
    prior,  # Prior information
    Samplesizemax = max_sample_resource,  # Maximum sample size
    pre_defined_ci = input$pre_defined_ci[i],  # Pre-defined credible interval
    groups_of_data = groups_of_data,  # Number of data groups
    MaxDatainfuture = max_sample_resource,  # Maximum data in future
    tolerance_level = input$tolerance_level[i],  # Tolerance level
    Nmin = 3,  # Minimum sample size for regression model
    nSim = 300  # Number of simulations
  )
  
  print(paste("Loop", i, "is finished"))  # Print loop completion message
  result[[i]] <- res  # Save result
}

# Save results for plotting later
save(result, file = paste("Results/", prior_type, ".RData", sep = ""))
