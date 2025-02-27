source("Rscripts/BOS.R")
library("readxl")
library("gam")
library("progress")

load("Data/central_informative_prior_statisticstandard.RData")
groups_of_data = 400
# group_n = 5
input <- expand.grid(pre_defined_ci = as.numeric(pre_defined_ci),
                     tolerance_level = seq(0, 1,0.1),
                     Nmin_for_regressionmodel = seq(10, 30, 10))

result= replicate(length(input$pre_defined_ci), data.frame())
index = 0
# # Create the directory string
# output_directory <- paste("result", groupname,
#                           regressionname,datafile,"_group_n",group_n, sep = "_")
# # Create the directory if it doesn't exist
# dir.create(output_directory, showWarnings = FALSE)
# Set up progress bar.
pb <- progress_bar$new(format = "[:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       result = length(input$pre_defined_ci),
                       clear = FALSE)      

for(i in 1:length(input$pre_defined_ci)){
    gc()
    index = index +1
    pb$tick()
    #file_path <- file.path(output_directory, paste0("Figure", index, ".jpeg"))
    res <- BOS(INDATA,
              experiment_data_groups,
              prior,
                             # Regression = regression_setting,
                             # regression_model_choice =1,
              Samplesizemax = max_sample_resource,
              pre_defined_ci = input$pre_defined_ci[i],
              groups_of_data = groups_of_data,
              MaxDatainfuture = max_sample_resource,
              tolerance_level = input$tolerance_level[i],
              Nmin = 3,
              Nmin_for_regressionmodel = input$Nmin_for_regressionmodel[i],
              nSim = 300,
              n_picks = 1e5,
              index = index,
              file_path)
    print(paste("Loop",i,"is finished"))
    result[[i]] <- res
}
# Save results for plotting later
save(result, file = paste("Results/",prior_type, ".RData", sep = ""))
