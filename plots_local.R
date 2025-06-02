setwd("C:/Users/alexander/Dropbox/MARIA-REPULSIVEWEIGHTS/paper")
library(mcclust.ext)
library(mcclust)
library(ggalt)
library(dplyr)
library(ggsci)
library(salso)
library(viridisLite)
library(tidyr)

all_alphas <- c( 1)
all_gammas <- c(0, 0.25, 0.5,  1, 98, 100, 101)
all_zetas <- c(0.001,  1,  3, 5, 98,99)
repulsive_grid <- expand.grid(all_alphas, all_gammas, all_zetas)
colnames(repulsive_grid) <- c("alpha", "gamma", "zeta")
n_save = 5e3;  n_burn = 5e3; n_thin = 1
# prior:
repulsive_grid_prior <- data.frame(repulsive_grid) %>%  mutate(row = row_number()) %>% 
  filter((gamma %in% c(100))) %>% 
  filter(!zeta%in% c(98, 3)  )

psm_prior_plot_order <-  repulsive_grid_prior[order(repulsive_grid_prior[, 2],repulsive_grid_prior[, 3],repulsive_grid_prior[, 1]), ] %>% 
  pull(row)


all_M_a.w <- data.frame()

for(run in psm_prior_plot_order){
  print(run)
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )

  # M_a  
  M_a.w <- cbind(repulsive_grid[run,],
                 M_a = readRDS(paste0("simstudy/fullrep/M_a/M_a_",filename,".rds")))%>%
    mutate(iter = 1:n_save)
  
  all_M_a.w <- rbind(all_M_a.w, M_a.w) 
  
}


all_repara.w <- data.frame()
for(run in psm_prior_plot_order){
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  gamma <- 
  zetas <- readRDS(paste0("simstudy/fullrep/zetas/zetas_",filename,".rds"))
 
  repara.w <- cbind(repulsive_grid[run,],
              gamma_out = readRDS(paste0("simstudy/fullrep/gammas/gammas_",filename,".rds")),
              zeta = readRDS(paste0("simstudy/fullrep/zetas/zetas_",filename,".rds")))%>%
               mutate(iter = 1:n_save)

    
  all_repara.w <- rbind(all_repara.w, repara.w) 
  
  
}

all_repara.l <- all_repara.w %>% 
                  pivot_longer(cols =-c("alpha", "gamma", "zeta", "iter"), values_to = "value", names_to = "repara")












# Plots
library(patchwork)
path <- "figures/fullrep/M_random/ratio/"



create_paper_plots <- function(df, geom, ncol, legend){
  
  
  if(geom == "line"){
   
    plot <-   df %>% 
      filter(repara == "gamma_out") %>% 
      mutate(repara = ifelse(repara == "gamma_out", "", repara)) %>% 
      filter(zeta != 3) %>% 
      mutate(zeta = ifelse(zeta == 99, 100, zeta)) %>% 
      ggplot(aes(x = iter, y = value, col = factor(zeta))) +
      geom_line()+
      facet_wrap(gamma ~ zeta, scales = "free", ncol = ncol) + 
      labs(title = "",
           x = "",
           y = "",
           col = expression(gamma)) +
      theme_minimal() +
      scale_color_npg() +
      if(legend == TRUE){
        theme(legend.position = "bottom") 
      }else{
        theme(legend.position = "none") 
      }  +
      theme(panel.spacing = unit(1, "cm"),  # Increase the distance between plots
            text=element_text(size=20), #change font size of all text
            strip.text.x = element_blank(),
            axis.text=element_text(size=15), #change font size of axis text
            axis.title=element_text(size=15), #change font size of axis titles
            plot.title=element_text(size=15), #change font size of plot title
            legend.text=element_text(size=15), #change font size of legend text
            legend.title=element_text(size=15)) 
    
  }else if(geom == "density"){
    
    plot <- df %>% 
      filter(repara == "gamma_out") %>% 
      filter(zeta != 3) %>% 
      mutate(zeta = ifelse(zeta == 99, 100, zeta)) %>% 
      ggplot(aes(x = value, fill = factor(zeta), col = factor(zeta))) +
      geom_density(alpha = 0.4)+
      facet_wrap(zeta ~ ., scales = "free", ncol = ncol) + 
      labs(title = "",
           x = "",
           y = "",
           col = expression(gamma),
           fill = expression(gamma)) +
      theme_minimal() +
      scale_fill_npg() +
      scale_color_npg() +
      if(legend == TRUE){
        theme(legend.position = "bottom") 
      }else{
        theme(legend.position = "none") 
      }  +
      theme(panel.spacing = unit(1, "cm"),  # Increase the distance between plots
            text=element_text(size=20), #change font size of all text
            #strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text=element_text(size=15), #change font size of axis text
            axis.title=element_text(size=15), #change font size of axis titles
            plot.title=element_text(size=15), #change font size of plot title
            legend.text=element_text(size=15), #change font size of legend text
            legend.title=element_text(size=15)) 
    
    
  }else if(geom == "barplot"){
    
  plot <-  df %>% 
      filter(zeta != 3) %>% 
      mutate(zeta = ifelse(zeta == 99, 100, zeta)) %>% 
      ggplot(aes(x = M_a, fill = factor(zeta)))+
      geom_bar(aes(y = after_stat(prop)),position = "dodge2") +
      facet_wrap(zeta ~., scales = "fixed", ncol = ncol) + 
      labs(title = "",
           x = "",
           y = "",
           fill = expression(zeta /gamma)) +
      theme_minimal() +
      scale_fill_npg() +
    # if(legend == TRUE){
    #   theme(legend.position = "bottom") 
    # }else{
    #   theme(legend.position = "none") 
    # }  +
      theme(panel.spacing = unit(1, "cm"), 
            text=element_text(size=20), #change font size of all text
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.text=element_text(size=15), #change font size of axis text
            axis.title=element_text(size=15), #change font size of axis titles
            plot.title=element_text(size=15), #change font size of plot title
            legend.text=element_text(size=15), #change font size of legend text
            legend.title=element_text(size=15))+
      scale_x_continuous(breaks = scales::breaks_pretty(n = 10)) 
    
  }

  return(plot)
  
}

# all_M_a.w
# all_repara.l


# M_bar

pdf(paste0(path,"M_a_bar_ratio_ncol4.pdf"), width = 16, height = 8)

plot <- create_paper_plots(all_M_a.w, "barplot", 4,legend = TRUE)

print(plot)

dev.off()

pdf(paste0(path,"M_a_bar_ratio_ncol1.pdf"), width = 8, height = 16)

plot <- create_paper_plots(all_M_a.w, "barplot", 1,legend = TRUE)

print(plot)

dev.off()


# geom_line

pdf(paste0(path,"gamma_line_ratio_ncol4.pdf"), width = 16, height = 8)

plot <- create_paper_plots(all_repara.l, "line", 4,legend = TRUE)

print(plot)

dev.off()

pdf(paste0(path,"gamma_line_ratio_ncol1.pdf"), width = 8, height = 16)

plot <- create_paper_plots(all_repara.l, "line", 1,legend = TRUE)

print(plot)

dev.off()

# geom_density

pdf(paste0(path,"gamma_density_ratio_ncol4.pdf"), width = 16, height = 8)

plot <- create_paper_plots(all_repara.l, "density", 4,legend = TRUE)

print(plot)

dev.off()

pdf(paste0(path,"gamma_density_ratio_ncol1.pdf"), width = 8, height = 16)

plot <- create_paper_plots(all_repara.l, "density", 1,legend = TRUE)

print(plot)

dev.off()

# Ncol 1 combined:

# Example ggplot objects
p1 <- create_paper_plots(all_repara.l, "line", 1,legend = FALSE) #+ theme(legend.position = "none")
p2 <-  create_paper_plots(all_repara.l, "density", 1,legend = FALSE) #+ theme(legend.position = "none")
p3 <- create_paper_plots(all_M_a.w, "barplot", 1,legend = TRUE) #+ theme(legend.position = "bottom")

pdf(paste0(path,"M_bar_repara_line_dens.pdf"), width = 24, height = 16)
# Combine plots with a shared legend
(p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + (p1 + theme(legend.position = "none"))   + (p2 + theme(legend.position = "none"))   
dev.off()

# Legend on the side:

pdf(paste0(path,"M_bar_repara_line_dens_legendside.pdf"), width = 24, height = 16)
# Combine plots with a shared legend
 (p1 + theme(legend.position = "none"))   + (p2 + theme(legend.position = "none")) +  (p3 + plot_layout(guides = "collect")) 
dev.off()

# PSM:

png( paste0(path, "psm_ratio_ncol4.png"),
     width     = 8,
     height    = 2,
     units     = "in",
     res       = 500,
     pointsize = 4)

layout(matrix(c(1:nrow(repulsive_grid_prior)), length(unique(repulsive_grid_prior$gamma)), length(unique(repulsive_grid_prior$zeta)),byrow = TRUE))
par(mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
for(run in psm_prior_plot_order){
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/fullrep/allocs/allocs_",filename,".rds"))
  #plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")
  psm <- psm(allocs)
  hc = hclust(as.dist(1 - psm), method = "complete", members = NULL)
  psm_hc = psm
  n = nrow(psm)
  psm_hc[1:n, ] = psm_hc[hc$order, ]
  psm_hc[, 1:n] = psm_hc[, hc$order]
  image(1:n, 1:n, 1 - psm_hc, col = viridis(20), xlab = "", ylab = "")
  
  
}

dev.off()




png( paste0(path,"psm_ratio_ncol1.png"),
     width     = 2,
     height    = 8,
     units     = "in",
     res       = 500,
     pointsize = 4)

layout(matrix(c(1:nrow(repulsive_grid_prior)), length(unique(repulsive_grid_prior$gamma)) * length(unique(repulsive_grid_prior$zeta)),1 ,byrow = TRUE))
par(mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
for(run in psm_prior_plot_order){
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/fullrep/allocs/allocs_",filename,".rds"))
  #plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")
  psm <- psm(allocs)
  hc = hclust(as.dist(1 - psm), method = "complete", members = NULL)
  psm_hc = psm
  n = nrow(psm)
  psm_hc[1:n, ] = psm_hc[hc$order, ]
  psm_hc[, 1:n] = psm_hc[, hc$order]
  image(1:n, 1:n, 1 - psm_hc, col = viridis(20), xlab = "", ylab = "")
  
}

dev.off()

png( paste0(path,"psm_ratio_ncol1_2.png"),
     width     = 4,
     height    = 16,
     units     = "in",
     res       = 500,
     pointsize = 4)

layout(matrix(c(1:nrow(repulsive_grid_prior)), length(unique(repulsive_grid_prior$gamma)) * length(unique(repulsive_grid_prior$zeta)),1 ,byrow = TRUE))
par(mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
for(run in psm_prior_plot_order){
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/fullrep/allocs/allocs_",filename,".rds"))
  #plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")
  psm <- psm(allocs)
  hc = hclust(as.dist(1 - psm), method = "complete", members = NULL)
  psm_hc = psm
  n = nrow(psm)
  psm_hc[1:n, ] = psm_hc[hc$order, ]
  psm_hc[, 1:n] = psm_hc[, hc$order]
  image(1:n, 1:n, 1 - psm_hc, col = viridis(20), xlab = "", ylab = "")
  
}

dev.off()

# bigger margins:

png( paste0(path,"psm_ratio_ncol1_largemarge.png"),
     width     = 4,
     height    = 16,
     units     = "in",
     res       = 500,
     pointsize = 4)

layout(matrix(c(1:nrow(repulsive_grid_prior)), length(unique(repulsive_grid_prior$gamma)) * length(unique(repulsive_grid_prior$zeta)),1 ,byrow = TRUE))
par(mar = c(2, 2, 5, 2), oma = c(0.5, 0.5, 0.5, 0.5))
for(run in psm_prior_plot_order){
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/fullrep/allocs/allocs_",filename,".rds"))
  #plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")
  psm <- psm(allocs)
  hc = hclust(as.dist(1 - psm), method = "complete", members = NULL)
  psm_hc = psm
  n = nrow(psm)
  psm_hc[1:n, ] = psm_hc[hc$order, ]
  psm_hc[, 1:n] = psm_hc[, hc$order]
  image(1:n, 1:n, 1 - psm_hc, col = viridis(20), xlab = "", ylab = "")
  
}

dev.off()



# Same as pdf:

pdf( paste0(path,"psm_ratio_ncol1_2.png"), width = 8, height = 16)

layout(matrix(c(1:nrow(repulsive_grid_prior)), length(unique(repulsive_grid_prior$gamma)) * length(unique(repulsive_grid_prior$zeta)),1 ,byrow = TRUE))
par(mar = c(2, 2, 2, 2), oma = c(0, 0, 0, 0))
for(run in psm_prior_plot_order){
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/fullrep/allocs/allocs_",filename,".rds"))
  #plotpsm(psm(allocs), main = paste0("PSM for ", gsub("_", " = ", filename)), method = "complete")
  psm <- psm(allocs)
  hc = hclust(as.dist(1 - psm), method = "complete", members = NULL)
  psm_hc = psm
  n = nrow(psm)
  psm_hc[1:n, ] = psm_hc[hc$order, ]
  psm_hc[, 1:n] = psm_hc[, hc$order]
  image(1:n, 1:n, 1 - psm_hc, col = viridis(20), xlab = "", ylab = "")
  
}

dev.off()



########################################################################
# PSMs done with ggplot
########################################################################


#pdf( paste0(path,"psm_ratio_ncol1_2.png"), width = 8, height = 16)
gg_psms <- list()
for(run in psm_prior_plot_order){
  
  filename <- paste0( paste0(colnames(repulsive_grid[run,]), "_"), paste0(repulsive_grid[run,]), collapse = "_" )
  allocs <- readRDS(paste0("simstudy/fullrep/allocs/allocs_",filename,".rds"))
  
ind <- which(run == psm_prior_plot_order)

# Compute the Posterior Similarity Matrix from the allocation matrix
posterior_similarity <- matrix(0, nrow = ncol(allocs), ncol = ncol(allocs))

for (i in 1:nrow(allocs)) {
  posterior_similarity <- posterior_similarity + (outer(allocs[i, ], allocs[i, ], "=="))
}

# Normalize the similarity matrix by the number of iterations
posterior_similarity <- posterior_similarity / nrow(allocs)

# Convert the matrix to a format suitable for ggplot2
library(reshape2)
similarity_df <- melt(posterior_similarity)
colnames(similarity_df) <- c("Observation1", "Observation2", "Similarity")

# Perform hierarchical clustering on the similarity matrix
hc <- hclust(as.dist(1 - posterior_similarity))
reordered_indices <- hc$order
similarity_df$Observation1 <- factor(similarity_df$Observation1, levels = reordered_indices)
similarity_df$Observation2 <- factor(similarity_df$Observation2, levels = reordered_indices)

# Very Light Gray: #e0e0e0
# Ivory: #fffff0
# Lavender: #f4f3ff
# Beige: #f5f5dc
# Mint Cream: #f5fff5


# good: "#f0f0f0"

# Plot the clustered heatmap
gg_psms[[ind]] <- ggplot(similarity_df, aes(x = Observation1, y = Observation2, fill = Similarity)) +
  geom_tile() +
  scale_fill_gradient(low = "#f4f3ff", high = pal_npg()(4)[ind]) +
  #scale_fill_gradient(low = "white", high = "darkred") +
  labs(title = "",
       x = "",
       y = "",
       fill = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

}

gg_psm_nps <- gg_psms[[1]] / gg_psms[[2]] / gg_psms[[3]] / gg_psms[[4]] 

# Ncol 1 combined now with psm
  
# Example ggplot objects
p1 <- create_paper_plots(all_repara.l, "line", 1,legend = FALSE) #+ theme(legend.position = "none")
p2 <-  create_paper_plots(all_repara.l, "density", 1,legend = FALSE) #+ theme(legend.position = "none")
p3 <- create_paper_plots(all_M_a.w, "barplot", 1,legend = TRUE) #+ theme(legend.position = "bottom")

# Legend on the side:

pdf(paste0(path,"M_bar_repara_line_dens_psm_legendside.pdf"), width = 32, height = 16)

png( paste0(path,"M_bar_repara_line_dens_psm_legendside_2.png"),
     width     = 32,
     height    = 16,
     units     = "in",
     res       = 500,
     pointsize = 4)

# Combine plots with a shared legend

(p1 + theme(legend.position = "none"))  | (p2 + theme(legend.position = "none"))| (gg_psm_nps + theme(legend.position = "none")) | (p3 + plot_layout(guides = "collect")) 

dev.off()

