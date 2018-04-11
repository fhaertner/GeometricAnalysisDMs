
library(dplyr)
library(purrr)
library(igraph)
library(NetHypGeom)


#' z-test
#' 
#'  @param x The value to compare against a random distribution
#'  @param rnd_distr The random distribution
#'  @param alt Whether to test for 'bigger' or 'smaller' than random
#'  
z_test <- function(x, rnd_distr, alt = "bigger"){
  z_score <- (x - mean(rnd_distr))/sd(rnd_distr)
  if(alt == "bigger"){
    return(pnorm(-z_score))
  }else{
    return(pnorm(z_score))
  }
}

#' Returns a single-row tibble with all stats about connectivity patterns of 
#' proteins associated with single disease in the hPIN
#' 
#' @param dm A tibble with genes associated with a single disease
#' @param net An igraph object representing the hPIN
#' @param H2 Pairwise hyperbolic distance matrix
#' @param epochs Number of random samplings
#' 
compute_dm_stats <- function(dm, net, H2, epochs = 1000){
  stats <- tibble(disease = dm$diseaseName[1], # Disease name
                  assoc = nrow(dm),            # Num. of associated genes
                  dm_size = integer(1),        # Num. of genes in the hPIN
                  comp = integer(1),           # Number of connected components
                  lcc = integer(1),            # Size of the larget component
                  pval_lcc = numeric(1),       # Is LCC larger than expected?
                  sp = numeric(1),             # Avg. SP to closest protein
                  pval_sp = numeric(1),        # Are genes closer than expected?
                  h2 = numeric(1),             # Avg. hyp. dist. to closest gene
                  pval_h2 = numeric(1))        # Are genes closer than expected 
                                               # in hyperbolic space?
  
  genes <- as.character(unique(dm$geneId))
  
  stats$assoc[1] <- length(genes)
  stats$dm_size[1] <- sum(genes %in% V(net)$name)
  
  # Only consider disease genes present in the hPIN
  genes <- V(net)$name[V(net)$name %in% genes]
  idx <- which(V(net)$name %in% genes)
  
  # Focus on the subgraph formed by the disease genes
  sub_net <- induced_subgraph(net, genes)
  comps <- components(sub_net)
  
  stats$comp[1] <- comps$no
  stats$lcc[1] <- max(comps$csize)
  
  sp_matrix <- distances(net, genes, genes)
  diag(sp_matrix) <- Inf
  stats$sp[1] <- mean(apply(sp_matrix, 1, min))
  
  h2_matrix <- H2[idx, idx]
  diag(h2_matrix) <- Inf
  stats$h2[1] <- mean(apply(h2_matrix, 1, min))
  
  # Comparisons against random sets of proteins
  rnd_lcc <- integer(epochs)
  rnd_sp <- integer(epochs)
  rnd_h2 <- numeric(epochs)
  
  # Form a pool of proteins that are not related to the disease but have
  # similar degrees
  # deg <- degree(net, genes)
  # pool <- which((degree(net) >= min(deg)) & (degree(net) <= max(deg)))
  # pool <- V(net)$name[setdiff(pool, idx)]
  
  for(i in 1:epochs){
    # Sample the number of disease genes at random from the hPIN
    rnd_genes <- sample(V(net)$name, length(genes))
    # rnd_genes <- sample(pool, length(genes))
    idx <- which(V(net)$name %in% rnd_genes)
    
    # Focus on the subgraph formed by these genes
    rnd_subn <- induced_subgraph(net, rnd_genes)
    
    rnd_lcc[i] <- max(components(rnd_subn)$csize)
    
    sp_matrix <- distances(net, rnd_genes, rnd_genes)
    diag(sp_matrix) <- Inf
    rnd_sp[i] <- mean(apply(sp_matrix, 1, min))
    
    h2_matrix <- H2[idx, idx]
    diag(h2_matrix) <- Inf
    rnd_h2[i] <- mean(apply(h2_matrix, 1, min)) 
  }

  # Compute the fraction of times the real values are bigger or smaller than
  # the random ones
  # stats$pval_lcc[1] <- sum(rnd_lcc >= stats$lcc[1])/epochs
  # stats$pval_sp[1] <- sum(rnd_sp <= stats$sp[1])/epochs
  # stats$pval_h2[1] <- sum(rnd_h2 <= stats$h2[1])/epochs
  
  # Statistically test whether real values are bigger or smaller than
  # the random ones
  stats$pval_lcc[1] <- z_test(stats$lcc[1], rnd_lcc, "bigger")
  stats$pval_sp[1] <- z_test(stats$sp[1], rnd_sp, "smaller")
  stats$pval_h2[1] <- z_test(stats$h2[1], rnd_h2, "smaller")
  
  return(stats)
}

load("data/DisGeNETv5.RData")
load("data/hint_coords.RData")
hint_net <- readRDS("data/hint.rds")

coords <- coords_lh$polar
H <- sapply(seq(nrow(coords)), function(i) hyperbolic_dist(coords[i, ], coords))
save(H, file = "data/hint_H.RData")

dm_stats <- disgenet %>% 
  split(.$diseaseName) %>% 
  map_dfr(compute_dm_stats, hint_net, H, 1000)

save(dm_stats, file = "results/dm_stats_hint.RData")

