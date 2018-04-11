
library(igraph)
library(dplyr)
library(ggplot2)
library(cowplot)

# Given a network and its hyperbolic coordinates, plot the angular distribution 
# of its protein classes
plot.theta.distr <- function(net, coords, classes, bins){
  class.names <- colnames(classes[, -1])
  
  thetas <- c()
  types <- c()
  for(i in 1:length(class.names)){
    # Determine which network proteins are in class
    idx <- which(V(net)$symbol %in% classes$symbol[classes[, class.names[i]]])
    thetas <- c(thetas, coords$theta[idx])
    types <- c(types, rep(class.names[i], length(idx)))
  }
  df <- data.frame(theta = thetas, 
                   class = factor(types, levels = unique(types), ordered = T))
  p <- ggplot(df, aes(theta, ..density.., fill = class)) + 
    geom_histogram(position = "identity", bins = bins) + 
    scale_fill_manual(values = c("#ED5051", "#FE880F", "#45A939", "#FB9A99", 
                                 "#87C196", "#A6CEE3")) +
    labs(x = expression(theta), y = "Density") + theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(), 
          legend.justification = c(0,1), legend.position = c(0,1))
  return(p)
}

load("data/protein_classes.RData")

# Bio stats for hPIN ------------------------------------------------------

load("data/hPIN_150k.RData")
load("data/coords_hPIN_150k.RData")

hpin.tht <- plot.theta.distr(hPIN, coords, 
                             protein.classes[, -c(5, 9, 10, 11)], bins = 75)

coords$age <- V(hPIN)$age
coords <- coords[!is.na(coords$age),]
hpin.age <- ggplot(coords, aes(factor(age, levels = 6:1, ordered = T), r)) + 
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(12, max(coords$r))) + 
  labs(x = "Protein age", y = "Node radial coordinate") + theme_bw()

save(hpin.age, hpin.tht, file = "results/bio_eval_hpin.RData")

# Bio stats for HINT ------------------------------------------------------

net <- readRDS("data/hint.rds")
load("data/hint_coords.RData")
coords <- coords_lh$polar

hint.tht <- plot.theta.distr(net, coords, 
                             protein.classes[, -c(5, 9, 10, 11)], bins = 60)

ages <- tibble(id = V(hPIN)$name, age = V(hPIN)$age)
coords <- left_join(coords, ages, by = "id")

coords <- coords[!is.na(coords$age),]
hint.age <- ggplot(coords, aes(factor(age, levels = 6:1, ordered = T), r)) + 
  geom_boxplot(outlier.shape = NA) + 
  coord_cartesian(ylim = c(10, max(coords$r))) + 
  labs(x = "Protein age", y = "Node radial coordinate") + theme_bw()

save(hint.age, hint.tht, file = "results/bio_eval_hint.RData")
