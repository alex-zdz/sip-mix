# Load and source
library("rdist")
library("dplyr")
library("plotly")
library("ggtern")
library("viridis")
library("ggplot2")
library(dplyr)

################################################################################
# 1. Ternary dots pwds
################################################################################
# points.df <- rbind(
#                     data.frame(A = 0.255, B = 0.3, C = 0.37, col = 1),
#                     data.frame(A = 1 / 9, B = 0.3, C = 1 / 3, col = 2),
#                     data.frame(A = 1 / 30, B = 0.3, C = 1 / 3, col = 3)
#                     #data.frame(A = 1 / 2000, B = 1 / 3, C = 1 / 3, col = "green3")
# )
# 
# ternary_dots_1 <- 
# ggtern(data = points.df, mapping = aes(x = A, y = B, z = C, col = factor(col))) +
#   geom_point( size = 3, show.legend = FALSE) +
#   geom_crosshair_tern(show.legend = FALSE) +
#   theme_light() +
#   scale_T_continuous(limits=c(0.0,1.0),
#                      breaks=seq(0,1,by=0.1),
#                      labels=seq(0,1,by=0.1)) + 
#   scale_L_continuous(limits=c(0.0,1.0),
#                      breaks=seq(0,1,by=0.1),
#                      labels=seq(0,1,by=0.1)) +
#   scale_R_continuous(limits=c(0.0,1.0),
#                      breaks=seq(0,1,by=0.1),
#                      labels=seq(0,1,by=0.1)) 
# 
#   
# ggtern::ggsave(filename = "ternary_dots_1.png", plot = ternary_dots_1, dpi = 300)


# Blue
# 0.46, 0.6, 0.5
# blue <- c(0.46, 0.6, 0.5)
# c(rdist(blue))
# prod(c(rdist(blue)))
# #prod(c(rdist(blue)))
# # red
# # 0.33, 0.27, 0.4
# red <- c(0.33, 0.27, 0.4)
# c(rdist(red))
# prod(c(rdist(red)))
# # green
# # 0.41, 0.15, 0.46
# green <- c(0.41, 0.15, 0.46)
# c(rdist(green))
# prod(c(rdist(green)))

# #ternary_dots_1 <- 
#   ggtern(mapping = aes(x = A, y = B, z = C)) +
#   #geom_point(size = 2) +
#   # geom_point(data = data.frame(A = 1 / 6, B = 1 / 3, C = 1 / 3), col = "red", size = 4) +
#   # geom_point(data = data.frame(A = 1 / 3, B = 1 / 6, C = 1 / 3), col = "red", size = 4) +
#   # geom_point(data = data.frame(A = 1 / 3, B = 1 / 3, C = 1 / 6), col = "red", size = 4) +
#   geom_point(data = data.frame(A = 0.255, B = 0.3, C = 0.37), col = "red2", size = 3) +
#   #geom_point(data = data.frame(A = 0.3, B = 0.255, C = 0.3), col = "red2", size = 3) +
#   #geom_point(data = data.frame(A = 0.3, B = 0.3, C = 0.255), col = "red2", size = 3) +
#   geom_point(data = data.frame(A = 1 / 9, B = 1 / 3, C = 1 / 3), col = "blue3", size = 3) +
#   geom_point(data = data.frame(A = 1 / 3, B = 1 / 9, C = 1 / 3), col = "blue3", size = 3) +
#   geom_point(data = data.frame(A = 1 / 3, B = 1 / 3, C = 1 / 9), col = "blue3", size = 3) +
#   geom_point(data = data.frame(A = 1 / 200, B = 1 / 3, C = 1 / 3), col = "green3", size = 3) +
#   geom_point(data = data.frame(A = 1 / 3, B = 1 / 200, C = 1 / 3), col = "green3", size = 3) +
#   geom_point(data = data.frame(A = 1 / 3, B = 1 / 3, C = 1 / 200), col = "green3", size = 3) +
  # scale_T_continuous(limits=c(0.0,1.0),
  #                    breaks=seq(0,1,by=0.1),
  #                    labels=seq(0,1,by=0.1)) +
  # scale_L_continuous(limits=c(0.0,1.0),
  #                    breaks=seq(0,1,by=0.1),
  #                    labels=seq(0,1,by=0.1)) +
  # scale_R_continuous(limits=c(0.0,1.0),
  #                    breaks=seq(0,1,by=0.1),
  #                    labels=seq(0,1,by=0.1)) +
#   geom_Tisoprop(value=c(.5),colour='grey', size = 1, linetype = "dashed") +
# 
#   #geom_Lisoprop(value=c(.5),colour='grey', size = 1, linetype = "dashed") +
#   #geom_Risoprop(value=c(.5),colour='grey', size = 1, linetype = "dashed") +
#   theme_rgbw()


# png("ternary_dots_1.png", width = 5480, height = 5480, units = "px")
# 
# print(ternary_dots_1)
# 
# dev.off()
#
#
#ggtern::ggsave(filename = "ternary_dots_1.png", plot = ternary_dots_1, dpi = 300)





################################################################################
# 2. Sample-based density plots for the Selberg Dirichlet
################################################################################
library("bayesm")
library("boot")

# sDir_mh <- function(sample_num, alpha, beta, gamma, starting_value, tuning_sd){
#   
#   M <- length(starting_value)
#   all_x <- matrix(NA, sample_num, M)
#   x_curr <-  starting_value
#   
#   for(g in 1:sample_num){ 
#     
#     x_prop <- bayesm::rdirichlet(alpha = c(rep(alpha, M-1), beta))
#     
#     log_acrate <- 2 * gamma * (
#       apply(as.matrix(rdist(x_prop)), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  )) - apply(as.matrix(rdist(x_curr)), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  ))
#     )
#     
#     # excluding the last element of w from the pwds:
#     # log_acrate <- 2 * gamma * (
#     #   apply(as.matrix(rdist(head(x_prop, -1))), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  )) - apply(as.matrix(rdist(head(x_curr, -1))), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  ))
#     # )
#     
#     log_acrate <- min(log_acrate, 1)
#     
#     if(log(runif(1)) < log_acrate){
#       x_curr <-  x_prop
#     }
#     
#     all_x[g, ] <-  x_curr
#     
#   }
#   
#   return(all_x)
#   
# }

# Unequal alphas:

sDir_mh <- function(sample_num, alphas, gamma, starting_value, tuning_sd){
  
  M <- length(starting_value)
  all_x <- matrix(NA, sample_num, M)
  x_curr <-  starting_value
  
  for(g in 1:sample_num){ 
    
    x_prop <- bayesm::rdirichlet(alpha = alphas)
    
    log_acrate <- 2 * gamma * (
      apply(as.matrix(rdist(x_prop)), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  )) - apply(as.matrix(rdist(x_curr)), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  ))
    )
    
    # excluding the last element of w from the pwds:
    # log_acrate <- 2 * gamma * (
    #   apply(as.matrix(rdist(head(x_prop, -1))), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  )) - apply(as.matrix(rdist(head(x_curr, -1))), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  ))
    # )
    
    log_acrate <- min(log_acrate, 1)
    
    if(log(runif(1)) < log_acrate){
      x_curr <-  x_prop
    }
    
    all_x[g, ] <-  x_curr
    
  }
  
  return(all_x)
  
}


sDir_mh_new <- function(sample_num, alphas, gamma, starting_value, tuning_sd){
  
  M <- length(starting_value)
  all_x <- matrix(NA, sample_num, M)
  x_curr <-  starting_value
  
  for(g in 1:sample_num){ 
    
    x_prop <- bayesm::rdirichlet(alpha = alphas)
    
    # log_acrate <- 2 * gamma * (
    #   apply(as.matrix(rdist(x_prop[-1])), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  )) - apply(as.matrix(rdist(x_curr[-1])), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  ))
    # )
    
    # excluding the last element of w from the pwds:
    log_acrate <- 2 * gamma * (
      apply(as.matrix(rdist(head(x_prop, -1))), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  )) - apply(as.matrix(rdist(head(x_curr, -1))), 1, function(x) sum(log(abs( x ))[!is.infinite(log(abs( x )))]  ))
    )
    
    log_acrate <- min(log_acrate, 1)
    
    if(log(runif(1)) < log_acrate){
      x_curr <-  x_prop
    }
    
    all_x[g, ] <-  x_curr
    
  }
  
  return(all_x)
  
}

library(rdist)

M <- 3
all_alphas <- c(0.5, 1)
all_gammas <- c(0, 0.5, 1, 3)
dir.create(paste0("tern_sdir/"), showWarnings = FALSE, recursive = TRUE)
ind = 0
all_terns_equal_a <- list()
for(gamma in all_gammas){
  for(alpha in all_alphas){
    
    ind = ind + 1
    print(paste0(ind))
    # 
    # filename <- paste0("tern_sdir_a_",alpha,"_gamma_",gamma, ".png")
    
     
      starting_value <- c(bayesm::rdirichlet(alpha =  rep(alpha, M)))
      
      #burn-in included
      burn_in <- 100
      n_plot <- 100000
      sample_num  <- burn_in + n_plot
      tuning_sd = 1
      sdir_samples <- sDir_mh(sample_num, alphas =  c(3, 3, 3), gamma = 1, starting_value, tuning_sd)
      sdir_samples_new <- sDir_mh_new(sample_num, alphas =  c(3, 3, 3), gamma = 1, starting_value, tuning_sd)
      
      
      
      
      # Ternary Plot
      colnames(sdir_samples) <- LETTERS[1:3]
      colnames(sdir_samples_new) <- LETTERS[1:3]
      # Save plot
      #png(filename, width = 500, height = 500)
      
      # png( paste0("tern_sdir/", filename),
      #      width     = 10,
      #      height    = 10,
      #      units     = "in",
      #      res       = 1200,
      #      pointsize = 4)
      
      # tmp <- 
      all_terns_equal_a[[ind]] <- 
       ggtern(data = sdir_samples, aes(x = A, y = B, z = C)) +
         geom_hex_tern(binwidth=0.015, show.legend = FALSE) +
         #geom_point(size=0.25, show.legend = FALSE) +
         scale_fill_viridis_c() +
         geom_Tisoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
         geom_Lisoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
         geom_Risoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
         scale_T_continuous(limits=c(0.0,1.0),
                            breaks=seq(0,1,by=0.1),
                            labels=seq(0,1,by=0.1)) +
         scale_L_continuous(limits=c(0.0,1.0),
                            breaks=seq(0,1,by=0.1),
                            labels=seq(0,1,by=0.1)) +
         scale_R_continuous(limits=c(0.0,1.0),
                            breaks=seq(0,1,by=0.1),
                            labels=seq(0,1,by=0.1)) +
          theme(plot.title = element_text(size=20),
                axis.text=element_text(size=20)) +
          labs(
            #title = bquote((alpha~","~beta~","~gamma) == (.(round(alpha,2))~","~.(round(beta,2))~","~.(round(gamma,2)))),
            T = "",
            L = "",
            R = ""
          )
       #dev.off() 
      tmp <- 
        ggtern(data = sdir_samples_new, aes(x = A, y = B, z = C)) +
        geom_hex_tern(binwidth=0.015, show.legend = FALSE) +
        #geom_point(size=0.25, show.legend = FALSE) +
        scale_fill_viridis_c() +
        geom_Tisoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
        geom_Lisoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
        geom_Risoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
        scale_T_continuous(limits=c(0.0,1.0),
                           breaks=seq(0,1,by=0.1),
                           labels=seq(0,1,by=0.1)) +
        scale_L_continuous(limits=c(0.0,1.0),
                           breaks=seq(0,1,by=0.1),
                           labels=seq(0,1,by=0.1)) +
        scale_R_continuous(limits=c(0.0,1.0),
                           breaks=seq(0,1,by=0.1),
                           labels=seq(0,1,by=0.1)) +
        theme(plot.title = element_text(size=20),
              axis.text=element_text(size=20)) +
        labs(
          #title = bquote((alpha~","~beta~","~gamma) == (.(round(alpha,2))~","~.(round(beta,2))~","~.(round(gamma,2)))),
          T = "",
          L = "",
          R = ""
        ) 
      
      tmp| all_terns_equal_a[[ind]]
      
    }
}

all_terns_equal_a[[3]]

################################################################################
# Equal and Unequal alphas
################################################################################


################################################################################
#26.09.2024: plot all at once faceted to keep same legend
################################################################################

################################################################################
#27.05.2025: Use SDir_mh new which excludes the last element
################################################################################
library(tidyr)

burn_in <- 100
n_plot <- 100000
sample_num  <- burn_in + n_plot
tuning_sd = 1
all_sdir <- data.frame()

# Equal alphas plots
all_alphas <- matrix(NA, 2, 3)
all_alphas[1, ] <- rep(0.5, 3)
all_alphas[2, ] <- rep(1, 3)

for(a in 1:nrow(all_alphas)){
  print(paste0(a))
  alphas <-   all_alphas[a,]
  alpha <- unique(alphas)
  
  for(g in c(0.5, 1, 3)){

#sdir_samples <- sDir_mh(sample_num, alphas =  alphas, gamma = g, starting_value, tuning_sd)
    sdir_samples <- sDir_mh_new(sample_num, alphas =  alphas, gamma = g, starting_value, tuning_sd)
    colnames(sdir_samples) <- LETTERS[1:3]

all_sdir <- sdir_samples %>% data.frame() %>% 
  mutate(alpha = alpha, gamma = g) %>% 
  #pivot_longer(cols = -c(alpha, gamma), values_to = "value", names_to = "component") %>% 
  rbind(all_sdir)

  }
  
}

#aes(fill = ..density..),

all_sdir %>% 
filter(alpha == 1, gamma == 1) %>% 
ggtern( aes(x = A, y = B, z = C)) +
  geom_hex_tern(binwidth=0.015,  aes(fill = after_stat(density))) +
  scale_fill_viridis_c() +
  geom_Tisoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
  geom_Lisoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
  geom_Risoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed")


ggtern(data = all_sdir, aes(x = A, y = B, z = C)) +
  geom_hex_tern(aes(fill = ..density..), binwidth=0.015) +
  #geom_hex_tern(aes(fill = ..density..),binwidth=0.015,  aes(fill = after_stat(density))) +
  scale_fill_viridis_c() +
  geom_Tisoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
  geom_Lisoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
  geom_Risoprop(value=c(.5),colour='grey50', size = 2, linetype = "dashed") +
  facet_grid(alpha~gamma) + 
  #facet_wrap(alpha~gamma) + 
  scale_T_continuous(limits=c(0.0,1.0),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  scale_L_continuous(limits=c(0.0,1.0),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  scale_R_continuous(limits=c(0.0,1.0),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  theme(plot.title = element_text(size=20),
        axis.text=element_text(size=20)) +
  labs(
    #title = bquote((alpha~","~beta~","~gamma) == (.(round(alpha,2))~","~.(round(beta,2))~","~.(round(gamma,2)))),
    T = "",
    L = "",
    R = ""
  )










################################################################################
# 02.10.2024 Dirichlet and Selberg Dirichlet 
################################################################################

# Adapt "all_alphas" matrix 

M <- 3

#dir.create(paste0("tern_sdir/"), showWarnings = FALSE, recursive = TRUE)
# The figures are then moved to paper/figures/tern_sdir/unequal_alpha
setwd("C:/Users/alexander/Dropbox/MARIA-REPULSIVEWEIGHTS/paper/figures/tern_sdir")

# "Unequal" alphas plots
all_alphas <- matrix(NA, 6, 3)

all_alphas[1, ] <- c(2, 5, 2)
all_alphas[2, ] <- c(5, 2, 5)
all_alphas[3, ] <- c(5, 5, 2)
all_alphas[4, ] <- c(5, 3, 2)

# Equal alphas plots
#all_alphas <- matrix(NA, 2, 3)
all_alphas[5, ] <- rep(0.5, 3)
all_alphas[6, ] <- rep(1, 3)

for(a in 1:nrow(all_alphas)){
  
  alphas <-   all_alphas[a,]
  
  if(length(unique(alphas)) == 1){
    all_gammas <-  c(0.5, 1, 3)
  }else{
    all_gammas <-  c(1)
  }
    
    for(g in all_gammas){
      
      #gamma = 2
      #filename <- paste0("tern_sdir_legend_a_", paste0(alphas, collapse = "_"),"_gamma_",g, ".png")
      filename <- paste0("tern_sdir_a_", paste0(alphas, collapse = "_"),"_gamma_",g, ".png")
      starting_value <- c(bayesm::rdirichlet(alpha =  alphas))
      
      #burn-in included
      burn_in <- 100
      n_plot <- 100000
      sample_num  <- burn_in + n_plot
      tuning_sd = 1
      sdir_samples <- sDir_mh(sample_num, alphas =  alphas, gamma = g, starting_value, tuning_sd)
      colnames(sdir_samples) <- LETTERS[1:3]
      
      
      # Save plot
      png( paste0(filename),
           width     = 5,
           height    = 5,
           units     = "in",
           res       = 500,
           pointsize = 4)
      
      tmp <- 
        ggtern(data = sdir_samples, aes(x = A, y = B, z = C)) +
        geom_hex_tern(  aes(fill = after_stat(density)) ,binwidth=0.015, show.legend = FALSE) +
        #geom_point(size=0.25, show.legend = FALSE) +
        scale_fill_viridis_c() +
        geom_Tisoprop(value=c(.5),colour='grey50', size = 1, linetype = "dashed") +
        geom_Lisoprop(value=c(.5),colour='grey50', size = 1, linetype = "dashed") +
        geom_Risoprop(value=c(.5),colour='grey50', size = 1, linetype = "dashed") +
        scale_T_continuous(limits=c(0.0,1.0),
                           breaks=seq(0,1,by=0.1),
                           labels=seq(0,1,by=0.1)) +
        scale_L_continuous(limits=c(0.0,1.0),
                           breaks=seq(0,1,by=0.1),
                           labels=seq(0,1,by=0.1)) +
        scale_R_continuous(limits=c(0.0,1.0),
                           breaks=seq(0,1,by=0.1),
                           labels=seq(0,1,by=0.1)) +
        labs(
          #title = bquote((alpha~","~beta~","~gamma) == (.(round(alpha,2))~","~.(round(beta,2))~","~.(round(gamma,2)))),
          T = "",
          L = "",
          R = "",
          fill = ""
        ) +
        theme(legend.position = "bottom",
              legend.text=element_text(angle = 45))
      #+
      # theme(text=element_text(size=20), #change font size of all text
      #       strip.background = element_blank(),
      #       strip.text.x = element_blank(),
      #       axis.text=element_text(size=15), #change font size of axis text
      #       axis.title=element_text(size=15), #change font size of axis titles
      #       plot.title=element_text(size=15), #change font size of plot title
      #       legend.text=element_text(size=15), #change font size of legend text
      #       legend.title=element_text(size=15)) #change font size of legend title
      
      
      
      
      
      print(tmp)
      
      dev.off() 
      
      
    }
    
 # }
  

#for(a in 1:nrow(all_alphas)){
# Standard Dircihlet distribution:
#filename <- paste0("tern_dir_legend_a_",paste0(alphas, collapse = "_"), ".png")
filename <- paste0("tern_dir_a_",paste0(alphas, collapse = "_"), ".png")

png( paste0(filename),
     width     = 5,
     height    = 5,
     units     = "in",
     res       = 500,
     pointsize = 4)

if(a == 5){
  dir_samples <- brms::rdirichlet(sample_num, alpha = rep(0.65,3))
}else{
 
  dir_samples <- brms::rdirichlet(sample_num, alpha = alphas)
}

colnames(dir_samples) <- LETTERS[1:3]

# ggtern(data = dir_samples, aes(x = A, y = B, z = C)) +
#   geom_hex_tern(binwidth=0.015, show.legend = FALSE) +
#   scale_fill_viridis_c()

if(a == 6){
  color_end = 0.1
}else{
  color_end = 1
}

tmp <- ggtern(data = dir_samples, aes(x = A, y = B, z = C)) +
  geom_hex_tern( aes(fill = after_stat(density)), binwidth=0.015, show.legend = FALSE) +
  #geom_point(size=0.25, show.legend = FALSE) +
  scale_fill_viridis_c(end = color_end) +
  geom_Tisoprop(value=c(.5),colour='grey50', size = 1, linetype = "dashed") +
  geom_Lisoprop(value=c(.5),colour='grey50', size = 1, linetype = "dashed") +
  geom_Risoprop(value=c(.5),colour='grey50', size = 1, linetype = "dashed") +
  scale_T_continuous(limits=c(0.0,1.0),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  scale_L_continuous(limits=c(0.0,1.0),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  scale_R_continuous(limits=c(0.0,1.0),
                     breaks=seq(0,1,by=0.1),
                     labels=seq(0,1,by=0.1)) +
  labs(
    #title = bquote((alpha~","~beta~","~gamma) == (.(round(alpha,2))~","~.(round(beta,2))~","~.(round(gamma,2)))),
    T = "",
    L = "",
    R = "",
    fill = ""
  )+
  theme(legend.position = "bottom",
        legend.text=element_text(angle = 45))
#  theme(plot.title = element_text(size=20),
#        axis.text=element_text(size=20)) +

print(tmp)

dev.off()
    

}



###################################################################################################

# Quick Check:

alphas <- rep(0.5, 3)
alphas <- rep(1, 3)
starting_value <- c(bayesm::rdirichlet(alpha =  alphas))
burn_in <- 100
n_plot <- 100000
sample_num  <- burn_in + n_plot
tuning_sd = 1
sdir_samples <- sDir_mh(sample_num, alphas =  alphas, gamma = 0.1, starting_value, tuning_sd)
colnames(sdir_samples) <- LETTERS[1:3]

sdir_gg <- ggtern(data = sdir_samples, aes(x = A, y = B, z = C)) +
  geom_hex_tern(binwidth=0.025, show.legend = FALSE) +
  scale_fill_viridis_c()


dir_samples <- brms::rdirichlet(sample_num, alpha = alphas)
colnames(dir_samples) <- LETTERS[1:3]

dir_gg <- ggtern(data = dir_samples, aes(x = A, y = B, z = C)) +
  geom_hex_tern(binwidth=0.025, show.legend = FALSE) +
  scale_fill_viridis_c()



sdir_gg | dir_gg




################################################################################
# 2.1 Standard Dirichlet density plot
################################################################################
# library(ggtext)
# clrs <- c(
#   "#FFBE00",  # MCRN yellow
#   "#B92F0A",  # MCRN red
#   "#792A26",  # MCRN maroon
#   "#54191B",  # MCRN brown
#   "#242424",  # MCRN dark gray
#   "#2660ae"   # Blue from MCR flag
# )
# 
# # First triangle: random points
# withr::with_seed(1234, {
#   draws_3_7_2 <- brms::rdirichlet(n = 1e5, alpha = c(3, 7, 2)) |> 
#     data.frame() |> 
#     set_names(c("x", "y", "z"))
# })
# 
# tern1 <- draws_3_7_2 |> 
#   ggtern(aes(x = x, y = y, z = z)) +
#   geom_point(size = 0.2, alpha = 0.1, color = clrs[5]) +
#   scale_L_continuous(breaks = 0:5 / 5, labels = 0:5 / 5, name = "α<sub>1</sub>") +
#   scale_T_continuous(breaks = 0:5 / 5, labels = 0:5 / 5, name = "α<sub>2</sub>") +
#   scale_R_continuous(breaks = 0:5 / 5, labels = 0:5 / 5, name = "α<sub>3</sub>") +
#   theme(
#     tern.axis.title.L = element_markdown(face = "bold", color = clrs[1], size = rel(1.2)),
#     tern.axis.title.T = element_markdown(face = "bold", color = clrs[2], size = rel(1.2)),
#     tern.axis.title.R = element_markdown(face = "bold", color = clrs[4], size = rel(1.2))
#   )
# 
# # Second triangle: actual densities
# # Plotting the results from ddirichlet() is more difficult than using dbeta() or
# # dnorm() or other univariate distributions. We can't just use geom_function().
# # Instead, we need to generate a dataset of all possible combinations of the
# # three columns (x, y, and z here), keep only the rows where they sum to one,
# # and then find the probability density values for those rows with ddirichlet().
# # It's a complex process, but it works :shrug:
# 
# # Create a sequence of values for x, y, and z
# values <- seq(0, 1, by = 0.005)
# 
# # Generate all possible combinations of x, y, and z that sum to 1
# df <- expand.grid(x = values, y = values, z = values) |>
#   filter(x + y + z == 1) |> 
#   rowwise() |>
#   mutate(density = pmap_dbl(
#     list(x, y, z),
#     ~ brms::ddirichlet(
#       as.numeric(c(..1, ..2, ..3)),
#       alpha = c(3, 7, 2)
#     )
#   )) |> 
#   filter(!is.nan(density))
# 
# tern2 <- ggtern(data = df, aes(x = x, y = y, z = z)) +
#   geom_point(aes(color = density))
#   # Mess with the breaks so that the gradient is more visible (basically flatten
#   # super high values)
#   # scale_color_gradientn(
#   #   colors = clrs[5:1], 
#   #   values = scales::rescale(x = c(0, 1, 3, 8, 13), from = c(0, 13)),
#   #   guide = "none"
#   # ) +
# 
# # ggtern objects don't work with {patchwork} or gridExtra::grid.arrange() or
# # {cowplot} or any of the plot-combining packages, but {ggtern} comes with its
# # own version of grid.arrange(), so we can use that
# ggtern::grid.arrange(tern1, tern2, ncol = 2)



library(ggsimplex)
library(ggplot2)

library(brms)
data = rdirichlet(n = 100, alpha = c(1,2,3))
data = as.data.frame(data)
colnames(data) = c("pmp_1", "pmp_2", "pmp_3")

df_dirichlet = data.frame(true_model = 1)
df_dirichlet$Alpha = list(c(1, 2, 3))


ggplot() +
  coord_fixed(ratio=1, xlim=c(0, 1), ylim=c(0, 1))+
  theme_void() +
  geom_simplex_canvas() + 
  stat_simplex_density(data=df_dirichlet, fun = ddirichlet,
                       args = alist(Alpha=Alpha))
################################################################################
# 3. Sample-based ( or not) density plots for the Gaussian Ensemble distribution
################################################################################


G <- function(p, zeta){
  
  ind <- 0:(p-1)
  G <- zeta^(-p / 2 - zeta * p * (p - 1) / 4) * (2 * pi) ^ (p / 2) *
    prod( gamma(1 + (1 + ind) * (zeta / 2)) / gamma(1 + zeta / 2) )  
  
  return(G)
  
}

ge_dens <- function(x, z){
  
  p <- length(x)
  dens <- (1 / G(p = p, zeta = z)) * prod(exp(- (z / 2) * x^2 ))  * prod(c(rdist(x)))^z
  return(dens)
  
}


N = 100
xy <- expand.grid(y1 = seq(-3, 3, length.out=N), y2 = seq(-3, 3, length.out=N))


dir.create(paste0("ge_dens/"), showWarnings = FALSE, recursive = TRUE)

zetas <- c(0.001, 0.1, 0.5, 1, 2, 3, 5, 10)
zetas <- c(0.1,  1, 2, 5)
#all_z <- data.frame()
for(z in zetas){
  
  z_tmp <- apply(xy, 1, function(x) ge_dens(x = x, z = z))
  z_tmp <- data.frame(xy, z_tmp) %>% mutate(zeta = z)
  #all_z <- rbind(all_z, z_tmp)
  
  filename <- paste0("bw_ge_dens_zeta_",z,".png")
  
  tmp <- z_tmp %>%  
    ggplot() + aes(x = y1, y = y2, z = z_tmp) +
    geom_contour_filled() +
    xlab("") + 
    ylab("") +
    theme_bw() +
    theme(text=element_text(size=20), #change font size of all text
          axis.text=element_text(size=30), #change font size of axis text
          axis.title=element_text(size=20), #change font size of axis titles
          plot.title=element_text(size=20), #change font size of plot title
          legend.text=element_text(size=25), #change font size of legend text
          legend.title=element_text(size=25),#change font size of legend title
          legend.position = "none")
    
  
  png( paste0("ge_dens/",filename),
       width     = 10,
       height    = 10,
       units     = "in",
       res       = 1200,
       pointsize = 4)
  
  print(tmp)
  
  dev.off()
  
  #ggsave(filename = filename, plot = tmp, dpi = 300)

}

# all_z %>%  filter(zeta == 1) %>% 
#   ggplot() + aes(x=y1,y=y2,z=z_tmp) +
#   geom_contour_filled() +
#   ggtitle(paste0("Zeta = ", z))











