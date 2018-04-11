
# Given an igraph network, it finds its Temperature and maps it to H2
# Run as follows:
# RScript map_to_h2.R netFile.RData outfile_name

library(NetHypGeom)

# Parsing arguments passed from command line ------------------------------

args <- commandArgs(trailingOnly = TRUE)
net <- readRDS(args[1])


# Find the network's temperature based on PS networks at Temp = 0 ---------

epochs <- 10

clust_at_zero <- numeric(epochs)

N <- vcount(net)
avg_k <- mean(degree(net))
gma <- 2.115#fit_power_law(degree(net))$alpha

for(i in 1:epochs){
  ps_net <- ps_model(N = N, avg.k = avg_k, gma = gma, Temp = 0)$network
  clust_at_zero[i] <- transitivity(ps_net, "average")
}

slope <- (0 - 1)/mean(clust_at_zero)
real_clust <- transitivity(net, "average")
Temp <- slope * real_clust + 1
#Temp <- 0.81

# Map to hyperbolic space

coords_lh <- labne_hm(net = net, gma = gma, Temp = Temp, 
                      k.speedup = 10, w = 2*pi)

save(coords_lh, Temp, file = args[2])


