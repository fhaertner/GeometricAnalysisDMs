
library(igraph)

load("data/hPIN_150k_old.RData")
load("data/coords_hPIN_150k_old.RData")
load("~/Desktop/Projects/disease_mod_geom/data/hPIN_150k_H_old.RData")

# Identify duplicates
dups <- which(duplicated(V(hPIN)$name))
dup.names <- V(hPIN)$name[dups]

# The duplicates of minimum degree are to be deleted
min.deg <- numeric(length = length(dup.names))

for(i in 1:length(min.deg)){
  idx <- which(V(hPIN)$name == dup.names[i])
  degs <- degree(hPIN, v = idx)
  min.deg[i] <- idx[which.min(degs)]
}

coords <- coords[-min.deg,]
H <- H[-min.deg, -min.deg]

save(coords, file = "data/coords_hPIN_150k.RData")
save(H, file = "~/Desktop/Projects/disease_mod_geom/data/hPIN_150k_H.RData")

# Rebuild the network
proteins <- data.frame(name = V(hPIN)$name, symbol = V(hPIN)$symbol, is.tf = V(hPIN)$is.tf, stringsAsFactors = F)
proteins <- proteins[!duplicated(proteins[, 1:2]), ]
edg <- as_data_frame(hPIN, what = "edges")
edg <- edg[!duplicated(edg), ]

hPIN <- graph_from_data_frame(edg, directed = F, vertices = proteins)
hPIN <- simplify(hPIN)

save(hPIN, file = "data/hPIN_150k.RData")
