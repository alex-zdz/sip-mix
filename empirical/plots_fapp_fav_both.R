# 01.04.2025 script to create all dashboards and plots for the results from the run with the food-avoidance and food-approach dataset

#setwd("~/paper/2d/eating_behavior_scripts")

#Load libraries
library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggalt)

# replicate repulsive grid used to run the analysis in "fapp_fav_both.R"
# Create grid
all_gammas <- c(0, 0.1, 0.5, 1, 2, 5)
all_zetas <- c(0.001, 0.1, 0.5, 1, 2, 3 )
#all_zetas <- c(0.001, 0.1, 0.5, 1, 2, 3, 5)
all_datasets <- c("bmi_pp_fa", "bmi_pp_fav")

# 01.04.2025 trying the food-avoidance food-approach dataset:
all_datasets <- c("fapp_fav")
all_datasets <- c("bmi_pp_fa")

all_updates <- c("fixed")
all_alphas <- c( 1, 3, 5) # 0. already run
repulsive_grid <- expand.grid(all_gammas, all_zetas, all_datasets, all_updates, all_alphas)
colnames(repulsive_grid) <- c("gamma", "zeta", "dataset", "update", "alpha")
all_datasets <- unique(repulsive_grid$dataset)
all_data <- list()

dataset <- all_datasets[1]

y <- readRDS(paste0("empirical/data/",dataset,".RDS"))
colnames(y) <- c("V1", "V2")
all_data[[1]] <- y

filename <- paste0( paste0(colnames(repulsive_grid[35,-3]), "_"), paste0(repulsive_grid[35,-3]), collapse = "_" )

n_save <- length(readRDS(paste0("empirical/results_",dataset,"/M_a/M_a_", filename,".rds")))
# Pre-determined length of evaluated grid:
post_pred_length <- sqrt(length(readRDS(paste0("empirical/results_",dataset,"/post_dens/post_dens_", filename,".rds"))))

from_to <- apply(y, 2, range)
y_line1 <- seq(from_to[1,1], from_to[2,1], length.out = 2e1)
y_line2 <- seq(from_to[1,2], from_to[2,2], length.out = 2e1)
y_grid1 <- expand.grid(y_line1, y_line2)

# all_M <- matrix(NA, nrow(repulsive_grid), n_save)
# for(run in 1:nrow(repulsive_grid)){
#   
#   # Assign variable based on run-id
#   zeta_samp <- repulsive_grid$zeta[run]
#   gamma_samp <- repulsive_grid$gamma[run]
#   dataset <- repulsive_grid$dataset[run]
#   update <- repulsive_grid$update[run]
#   
#   filename <- paste0( paste0(colnames(repulsive_grid[run,-3]), "_"), paste0(repulsive_grid[run,-3]), collapse = "_" )
#   
#   all_M[run, ] <- readRDS(paste0("empirical/results_",dataset,"/M_a/M_a_", filename,".rds"))
#   
# }
# 
# M_grid <- cbind(repulsive_grid, all_M)
# M.l <- M_grid %>% pivot_longer(cols = -c(1:ncol(repulsive_grid)), names_to = "iter", values_to = "M_a") %>% mutate(iter = as.numeric(iter))



############################################################################################
# Results
############################################################################################

dir.create(paste0("empirical/results_",dataset,"/dashboards/"), showWarnings = FALSE, recursive = TRUE)

all_zetas <- unique(repulsive_grid$zeta)
pps <- list()
binder.ggs <- list()
vi.ggs <- list()
for(run in 1:nrow(repulsive_grid)){
  
  print(run)
  
  zeta_samp <- repulsive_grid$zeta[run]
  gamma_samp <- repulsive_grid$gamma[run]
  dataset <- repulsive_grid$dataset[run]
  update <- repulsive_grid$update[run]
  alpha_samp <- repulsive_grid$alpha[run]
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,-3]), "_"), paste0(repulsive_grid[run,-3]), collapse = "_" )
  
  # post_pred <-  readRDS(paste0("empirical/results_",dataset,"/post_dens/post_dens_", filename,".rds"))
  # 
  # df <- data.frame(x = y_grid1$Var1, y = y_grid1$Var2, density = post_pred)
  # 
  # pps[[run]] <- ggplot(df, aes(x = x, y = y)) +
  #   geom_tile(aes(fill = density)) +
  #   geom_point(data = as.data.frame(y), aes(x = V1, y = V2 ), color = "black") +
  #   geom_contour(aes(z = density, color = ..level..), bins = 10, size = 0.75) +
  #   scale_fill_viridis_c()+
  #   scale_color_viridis_c(option = "plasma") +
  #   labs(title = paste0("Posterior predictive for zeta = ", zeta_samp, " and gamma = ", gamma_samp, " alpha = ", alpha_samp))
  # 
  ############################################################################################
  # Binder and VI
  ############################################################################################
  allocs_samp = readRDS(paste0("empirical/results_",dataset,"/allocs/allocs_", filename,".rds")) 
  
  # VI
  
  # s_vi <- c(salso(allocs_samp, loss=salso::VI(a=1)))
  # 
  # 
  # df <- data.frame(x = y[,1], y = y[,2], vi = s_vi) 
  # 
  # vi.ggs[[run]] <- 
  #   ggplot() +
  #   geom_point(data = df, aes(x = x, y = y, color = factor(vi)),  size = 2) +
  #   geom_encircle(data = df, aes(x = x, y = y, group = factor(vi), fill = factor(vi)),
  #                 size = 0.5, expand = 0, alpha = 0.15, s_shape=0.5,  spread=0.025) +
  #   labs(title = paste0("VI for zeta = ", zeta_samp, " and gamma = ", gamma_samp, " alpha = ", alpha_samp), x = "", y = "") +
  #   theme_minimal() +
  #   theme(text=element_text(size=20), #change font size of all text
  #         axis.text=element_text(size=30), #change font size of axis text
  #         axis.title=element_text(size=20), #change font size of axis titles
  #         plot.title=element_text(size=20), #change font size of plot title
  #         legend.text=element_text(size=25), #change font size of legend text
  #         legend.title=element_text(size=25),#change font size of legend title
  #         legend.position = "none") 
  
  
  
  ##############################################################################
  # Save individual plots
  ##############################################################################
  
  # pdf(paste0("empirical/results_",dataset,"/dashboards/", "pps_gamma_",gamma_samp,"_zeta_",zeta_samp,".pdf"), width = 12, height = 10)
  # print(pps[[run]])
  # dev.off()
  # 
  # pdf(paste0("empirical/results_",dataset,"/dashboards/", "vi_gamma_",gamma_samp,"_zeta_",zeta_samp,".pdf"), width = 12, height = 10)
  # print(vi.ggs[[run]])
  # dev.off()
  
  ##############################################################################
  # End individual plots
  ##############################################################################
  
  # Binder
  s_binder <-  c(salso(allocs_samp, loss=salso::binder(a=1)))
  
  df <- data.frame(x = y[,1], y = y[,2], Cluster = as.factor(s_binder))
  
  # binder.ggs[[run]] <- 
  #   ggplot() +
  #   geom_point(data = df, aes(x = x, y = y, color = factor(color)),  size = 2) +
  #   geom_encircle(data = df, aes(x = x, y = y, group = color, fill = factor(color)),
  #                 size = 0.5, expand = 0, alpha = 0.15, s_shape=0.5,  spread=0.025) +
  #   labs(title = paste0("Binder estimate for zeta = ", zeta_samp, " and gamma = ", gamma_samp), x = "BMI", y = "Food approach person parameter") +
  #   theme_minimal() 
  
 binder.ggs[[run]] <- 
   ggplot() +
   geom_point(data = df, aes(x = x, y = y, color = Cluster),  size = 2) +
   geom_encircle(data = df, aes(x = x, y = y, group = Cluster, fill = Cluster),
                 size = 0.5, expand = 0, alpha = 0.15, s_shape=0.5,  spread=0.025) +
     scale_fill_npg() +
     scale_color_npg() +
   labs(title = paste0("Binder zeta = ", zeta_samp, " gamma = ", gamma_samp, " alpha = ", alpha_samp), x = "BMI", y = "Food approach person parameter", legend = "test") + #
   theme_minimal() +
   theme(text=element_text(size=10), #change font size of all text
         axis.text=element_text(size=10), #change font size of axis text
         axis.title=element_text(size=10), #change font size of axis titles
         plot.title=element_text(size=10), #change font size of plot title
         legend.text=element_text(size=5), #change font size of legend text
         legend.title=element_text(size=5),#change font size of legend title
   ) 
 
   
 pdf(paste0("empirical/results_",dataset,"/dashboards/", "binder_gamma_",gamma_samp,"_zeta_",zeta_samp,".pdf"), width = 12, height = 10)
  print(binder.ggs[[run]])
 dev.off()
  
  # M_bar
  
 # M_bar_gfacet <- M.l %>% 
 #   filter(dataset == !!dataset, alpha == 0.5) %>% 
 #   mutate(zeta = as.factor(zeta), gamma = as.factor(gamma), iter = as.numeric(iter), alpha = as.factor(alpha)) %>% 
 #   ggplot(aes(x =  M_a, fill =  zeta)) + 
 #   geom_bar(width = 0.5, aes(y = after_stat(prop)), position = position_dodge2(width = 0.5, preserve = "single")) + 
 #   scale_x_continuous( breaks = min(M.l$M_a):max(M.l$M_a)) +
 #   facet_wrap(gamma~.) +
 #   theme_bw() +
 #   ggtitle(paste0("")) +
 #   xlab("")+
 #   ylab("") +
 #   theme(text=element_text(size=20), #change font size of all text
 #         axis.text=element_text(size=20), #change font size of axis text
 #         axis.title=element_text(size=20), #change font size of axis titles
 #         plot.title=element_text(size=20), #change font size of plot title
 #         legend.text=element_text(size=25), #change font size of legend text
 #         legend.title=element_text(size=25)) #change font size of legend title
 # 
 # M_bar_zfacet <- M.l %>% 
 #   filter(dataset == !!dataset, alpha == 0.5) %>% 
 #   mutate(zeta = as.factor(zeta), gamma = as.factor(gamma), iter = as.numeric(iter), alpha = as.factor(alpha)) %>% 
 #   ggplot(aes(x =  M_a, fill =  gamma)) + 
 #   geom_bar(width = 0.5, aes(y = after_stat(prop)), position = position_dodge2(width = 0.5, preserve = "single")) + 
 #   scale_x_continuous( breaks = min(M.l$M_a):max(M.l$M_a)) +
 #   facet_wrap(zeta~.) +
 #   theme_bw() +
 #   ggtitle(paste0("")) +
 #   xlab("")+
 #   ylab("") +
 #   theme(text=element_text(size=20), #change font size of all text
 #         axis.text=element_text(size=20), #change font size of axis text
 #         axis.title=element_text(size=20), #change font size of axis titles
 #         plot.title=element_text(size=20), #change font size of plot title
 #         legend.text=element_text(size=25), #change font size of legend text
 #         legend.title=element_text(size=25)) #change font size of legend title
 
}

# pdf(paste0("empirical/results_",dataset,"/dashboards/", "M_bar_gfacet.pdf"), width = 12, height = 10)
# print(M_bar_gfacet)
# dev.off()
# 
# pdf(paste0("empirical/results_",dataset,"/dashboards/", "M_bar_zfacet.pdf"), width = 12, height = 10)
# print(M_bar_zfacet)
# dev.off()
# 

pdf(paste0("empirical/results_",dataset,"/dashboards/", "binder.pdf"), width = 12, height = 10)

for(start_index in seq(1, nrow(repulsive_grid), by = 6)){
  do.call(grid.arrange, c(c(binder.ggs[start_index:(start_index + 5)]), ncol=3))
}

dev.off()



# pdf(paste0("empirical/results_",dataset,"/dashboards/", "vi.pdf"), width = 12, height = 10)
# 
# for(start_index in seq(1, length(pps), by = 6)){
#   do.call(grid.arrange, c(c(vi.ggs[start_index:(start_index + 5)]), ncol=3))
# }
# 
# dev.off()


# too many clusters
# zeta = 0.5, gamma = 0
# too much smoothing
# zeta = 0.5, gamma = 2
# right balance:
# zeta = 2, gamma = 2
library(patchwork)
pdf(paste0("empirical/results_",dataset,"/dashboards/", "test.pdf"), width = 12, height = 6)

#print(binder.ggs[[85]] + binder.ggs[[89]] + binder.ggs[[101]])

print(binder.ggs[[85]] + binder.ggs[[89]] + binder.ggs[[101]])

dev.off()









