# 01.04.2025 script to create all dashboards and plots for the results from the run with the food-avoidance and food-approach dataset

#setwd("~/paper/2d/eating_behavior_scripts")

#Load libraries
library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggalt)
library(salso)
# replicate repulsive grid used to run the analysis in "fapp_fav_both.R"
# Create grid
all_gammas <- c(0, 0.5, 1, 2)
all_zetas <- c(0.001, 0.1, 0.5, 1, 2)
all_datasets <- c("m_bmi_ch_bmi_w07", "m_bmi_ch_bmi_w10", "m_bmi_ch_bmi_w12", "m_bmi_ch_bmi_w14")
all_updates <- c("fixed")
all_alphas <- c(1, 5)
repulsive_grid <- expand.grid(all_gammas, all_zetas, all_datasets, all_updates, all_alphas)
colnames(repulsive_grid) <- c("gamma", "zeta", "dataset", "update", "alpha")


# Plots

dir.create(paste0("empirical/results_",dataset,"/dashboards/"), showWarnings = FALSE, recursive = TRUE)

binder.ggs <- list()
vi.ggs <- list()
for(run in 1:nrow(repulsive_grid)){
  
  print(run)

  zeta_samp <- repulsive_grid$zeta[run]
  gamma_samp <- repulsive_grid$gamma[run]
  dataset <- repulsive_grid$dataset[run]
  update <- repulsive_grid$update[run]
  alpha_samp <- repulsive_grid$alpha[run]
  
  y <- readRDS(paste0("empirical/data/",dataset,".RDS"))
  colnames(y) <- c("x", "V2")

  filename <- paste0( paste0(colnames(repulsive_grid[run,-3]), "_"), paste0(repulsive_grid[run,-3]), collapse = "_" )
  
  ############################################################################################
  # Binder and VI
  ############################################################################################
  allocs_samp = readRDS(paste0("empirical/results_",dataset,"/allocs/allocs_", filename,".rds")) 
  
  # VI
  
  s_vi <- c(salso(allocs_samp, loss=salso::VI(a=1)))
  
  
  df <- data.frame(x = y[,1], y = y[,2], vi = s_vi) 
  
  vi.ggs[[run]] <- 
    ggplot() +
    geom_point(data = df, aes(x = x, y = y, color = factor(vi)),  size = 2) +
    geom_encircle(data = df, aes(x = x, y = y, group = factor(vi), fill = factor(vi)),
                  size = 0.5, expand = 0, alpha = 0.15, s_shape=0.5,  spread=0.025) +
    labs(title = paste0("VI for dataset = ", dataset, " zeta = ", zeta_samp, " gamma = ", gamma_samp, " alpha = ", alpha_samp), x = "", y = "") +
    theme_minimal() +
    theme(text=element_text(size=20), #change font size of all text
          axis.text=element_text(size=30), #change font size of axis text
          axis.title=element_text(size=20), #change font size of axis titles
          plot.title=element_text(size=20), #change font size of plot title
          legend.text=element_text(size=25), #change font size of legend text
          legend.title=element_text(size=25),#change font size of legend title
          legend.position = "none") 
  
  # Binder
  s_binder <-  c(salso(allocs_samp, loss=salso::binder(a=1)))
  
  df <- data.frame(x = y[,1], y = y[,2], Cluster = as.factor(s_binder))
  
  binder.ggs[[run]] <- 
    ggplot(df) +
    geom_point(aes(x = df[[1]], y = df[[2]], color = Cluster), size = 2) +
    geom_encircle(aes(x = df[[1]], y = df[[2]], group = Cluster, fill = Cluster),
                  size = 0.5, expand = 0, alpha = 0.15, s_shape=0.5, spread=0.025) +
    labs(title = paste0("Binder for dataset = ", dataset, " zeta = ", zeta_samp, " gamma = ", gamma_samp, " alpha = ", alpha_samp), x = "Food approach person parameter", y = "Food avoidance person parameter", legend = "test") +
    theme_minimal() +
    theme(text=element_text(size=10), #change font size of all text
          axis.text=element_text(size=10), #change font size of axis text
          axis.title=element_text(size=10), #change font size of axis titles
          plot.title=element_text(size=10), #change font size of plot title
          legend.text=element_text(size=5), #change font size of legend text
          legend.title=element_text(size=5) #change font size of legend title
    ) 
  
}


# for(dataset in all_datasets){

# pdf(paste0("empirical/results_",dataset,"/dashboards/", "binder.pdf"), width = 12, height = 10)

# for(start_index in seq(1, length(pps), by = 6)){
#   do.call(grid.arrange, c(c(binder.ggs[start_index:(start_index + 5)]), ncol=3))
# }

# dev.off()

# pdf(paste0("empirical/results_",dataset,"/dashboards/", "vi.pdf"), width = 12, height = 10)

# for(start_index in seq(1, length(pps), by = 6)){
#   do.call(grid.arrange, c(c(vi.ggs[start_index:(start_index + 5)]), ncol=3))
# }

# dev.off()

# }

# Calculate the number of plots per dataset
plots_per_dataset <- nrow(repulsive_grid) / length(all_datasets)

for(dataset_index in seq_along(all_datasets)){
  dataset <- all_datasets[dataset_index]
  
  # Calculate the start and end indices for the current dataset's plots
  start_index <- (dataset_index - 1) * plots_per_dataset + 1
  end_index <- dataset_index * plots_per_dataset
  
  # Create a PDF for each dataset
  pdf(paste0("empirical/results_", dataset, "/dashboards/", "binder_", dataset, ".pdf"), width = 12, height = 10)
  
  # Plot the binder.ggs for the current dataset
  for(plot_index in seq(start_index, end_index, by = 6)){
    do.call(grid.arrange, c(c(binder.ggs[plot_index:min(plot_index + 5, end_index)]), ncol=3))
  }
  
  dev.off()
}


# pdf(paste0("empirical/results_",dataset,"/dashboards/", "test.pdf"), width = 12, height = 6)

# print(binder.ggs[[85]] + binder.ggs[[89]] + binder.ggs[[101]])

# dev.off()

