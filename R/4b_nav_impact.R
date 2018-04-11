library(NetHypGeom)
library(dplyr)
library(purrr)

#' Given a target distribution of node degrees and a pool of degrees, returns
#' a sample of the pool with a distribution similar to the target one
#' 
#' @param trg The target node degree distribution
#' @param pool The pool of node degrees
#' 
get_nodes_with_sim_degs <- function(trg, pool){
   idx <- which(names(pool) %in% names(trg))
   pool[idx] <- NA
   trg_dens <- approxfun(trg)
   w <- trg_dens(pool)
   w[is.na(w)] <- 0
   sampling <- sample(pool, length(trg), replace = FALSE, prob = w)
   return(which(names(pool) %in% names(sampling)))
}

#' Given a tibble with genes associated with a disease, the coordinates of their
#' products in hyperbolic space and a protein network, study impact of faulty 
#' disease proteins on navigability
#' 
#' @param df The tibble of gene-disease associations
#' @param coords The coordinates of gene products in hyperbolic space
#' @param net An igraph object representing the human protein network
#' @param st Number of source-targets to consider
#' @param faulty Number of faulty proteins to consider
#' 
impact_on_nav <- function(df, coords, net, st = 1000, 
                          epochs = 100, faulty = 20){
   # disease_idx <- which(df$geneId %in% V(net)$name)
   disease_idx <- which(V(net)$name %in% df$geneId)
   res <- tibble(eff_dis = numeric(epochs), # efficiency of routeings if disease proteins are marked as faulty
                 hop_dis = numeric(epochs), # average hop distance of the routetings with faulty disease proteins
                 eff_rnd = numeric(epochs), # efficiency of routeings if random proteins are marked as faulty
                 hop_rnd = numeric(epochs), # average hop distance of the routetings with faulty random proteins
                 faulty_dis = numeric(epochs), # list containing the gene id's of the faulty disease proteins
                 faulty_rnd = numeric(epochs)) # list containing the gene id's of the faulty random proteins
   
   for(i in 1:epochs){
      rnd_idx <- get_nodes_with_sim_degs(degree(net, V(net)$name[disease_idx]),
                                         degree(net))
      # Determine the source and targets using linear indexing
      idx <- sample(vcount(net) * vcount(net), st)
      src <- ((idx - 1) %% vcount(net)) + 1
      trg <- floor((idx - 1) / vcount(net)) + 1
      
      # sample faulty proteins
      faulty_dis <- sample(disease_idx, faulty)
      faulty_rnd <- sample(rnd_idx, faulty)
      # Perform the routeing
      eff_disease <- greedy_route_packets(net, coords, src, trg, 
                                          faulty_dis)
      eff_rnd <- greedy_route_packets(net, coords, src, trg, 
                                      faulty_rnd)
      
      res$eff_dis[i] <-  sum(eff_disease > 0)/st
      res$hop_dis[i] <-  mean(eff_disease[eff_disease > 0])
      res$eff_rnd[i] <-  sum(eff_rnd > 0)/st
      res$hop_rnd[i] <- mean(eff_rnd[eff_rnd > 0])
      res$faulty_dis[i] <- list(V(net)$name[faulty_dis])
      res$faulty_rnd[i] <- list(V(net)$name[faulty_rnd])
   }
   
   return(res)
}

load("../data/hPIN_150k.RData")
load("../data/coords_hPIN_150k.RData")
load("../data/DisGeNETv5.RData")

# Measure the reference network efficiency
num_to_select <- 10000
idx <- sample(vcount(hPIN) * vcount(hPIN), num_to_select)
src <- ((idx - 1) %% vcount(hPIN)) + 1
trg <- floor((idx - 1) / vcount(hPIN)) + 1
pin_routeing <- greedy_route_packets(hPIN, coords, src, trg)
ref_eff <- sum(pin_routeing > 0)/num_to_select
ref_hop <- mean(pin_routeing[pin_routeing > 0 & !is.infinite(pin_routeing)])

# Measure efficiency if disease and random faulty proteins are introduced
eff <- disgenet %>% 
  split(.$diseaseName)

#idx <- sample(1:157, 1)
#eff <- eff[idx]

eff <- eff %>% 
  map(impact_on_nav, coords, hPIN, st = 500, epochs = 10, faulty = 20)

# Finally measure impact on efficiency
# imp <- tibble(imp_eff_dis = numeric(nrow(disgenet)), 
#               imp_hop_dis = numeric(nrow(disgenet)),
#               imp_eff_rnd = numeric(nrow(disgenet)),
#               imp_hop_rnd = numeric(nrow(disgenet)))
imp <- tibble(imp_eff_dis = numeric(length(eff)), 
              imp_hop_dis = numeric(length(eff)),
              imp_eff_rnd = numeric(length(eff)),
              imp_hop_rnd = numeric(length(eff)),
              pval_eff = numeric(length(eff)),
              pval_hop = numeric(length(eff)))

imp$imp_eff_dis = map_dbl(eff, function(x) mean(x$eff_dis) - ref_eff)
imp$imp_hop_dis = map_dbl(eff, function(x) ref_hop - mean(x$hop_dis))
imp$imp_eff_rnd = map_dbl(eff, function(x) mean(x$eff_rnd) - ref_eff)
imp$imp_hop_rnd = map_dbl(eff, function(x) ref_hop - mean(x$hop_rnd))
imp$pval_eff = map_dbl(eff, function(x) 
  wilcox.test(x$eff_dis, x$eff_rnd, alternative = "less")$p.value)
imp$pval_hop = map_dbl(eff, function(x) 
  wilcox.test(x$hop_dis, x$hop_rnd, alternative = "greater")$p.value)

save(ref_eff, ref_hop, eff, imp, file = "../results/nav_impact.RData")
   