library(NetHypGeom)
library(igraph)
library(ggplot2)

get_hdist_vs_go <- function(polar, go, bins = 20){
  
  # The GO semantic similarity matrix
  go[lower.tri(go, diag = T)] <- 0
  
  # Compute the matrix of pairwise hyperbolic distances between nodes
  # WARNING: This can be quite computationally intensive!!!
  dist <- sapply(seq(nrow(polar)), function(i) hyperbolic_dist(polar[i, ], polar))
  
  max.dist <- max(dist)
  dist[lower.tri(dist, diag = T)] <- Inf
  
  steps <- seq(0, max.dist, length.out = bins)
  
  res <- data.frame(dist = steps, 
                    avg_go = vector("numeric", length = bins), 
                    sd_go =  vector("numeric", length = bins), 
                    se_go =  vector("numeric", length = bins),
                    stringsAsFactors = F)
  
  for(i in 1:(bins - 1)){
    res$avg_go[i] <- mean(go[(dist >= steps[i]) & (dist < steps[i + 1])])
    res$sd_go[i] <- sd(go[(dist >= steps[i]) & (dist < steps[i + 1])])
    res$se_go[i] <- sd(go[(dist >= steps[i]) & (dist < steps[i + 1])])/sqrt(length(go[(dist >= steps[i]) & (dist < steps[i + 1])]))
  }
  
  # Compute the last GO Sem Sims
  dist[lower.tri(dist, diag = T)] <- 0
  res$avg_go[bins] <- mean(go[dist >= steps[bins]])
  res$sd_go[bins] <- sd(go[dist >= steps[bins]])
  res$se_go[bins] <- sd(go[dist >= steps[bins]])/sqrt(length(go[dist >= steps[bins]]))
  
  return(res)
}

load("data/hPIN_150k.RData")
load("data/coords_hPIN_150k.RData")
load("data/bpSim_IEA.RData")
go <- bp.sim
rm(bp.sim)
# load("data/ccSim_IEA.RData")
# go <- cc.sim
# rm(cc.sim)

# Reduce the coords to the proteins present in the GO sim matrix
idx <- V(hPIN)$symbol %in% row.names(go)
coords <- coords[idx, ]

# Shrink the GO sim matrix
idx <- row.names(go) %in% V(hPIN)$symbol
go[idx, idx]

# Re-arrange the GO sim matrix to coincide with network order
idx <- V(hPIN)$symbol %in% row.names(go)
go <- go[V(hPIN)$symbol[idx], V(hPIN)$symbol[idx]]

# Perform the analysis
res <- get_hdist_vs_go(coords, go, 21)
p <- ggplot(res, aes(dist, avg_go)) + 
  geom_pointrange(aes(ymin = avg_go - se_go, ymax = avg_go + se_go)) +
  labs(x = "Hyperbolic distance", y = "GO semantic similarity") +
  theme_bw()

save(res, p, file = "results/h2_vs_bp.RData")
# save(res, p, file = "results/h2_vs_cc.RData")