
# Linear and non-linear dimensionality reduction applied to the 
# disease separation matrix

library(igraph)
library(dplyr)
library(ggplot2)
library(stringr)
library(RSpectra)

ncmce <- function(D, d = 2){
  
  # Form a fully-connected graph
  g <- graph_from_adjacency_matrix(D, "undirected", weighted = TRUE)
  
  # Compute the network's minimum spanning tree (MST)
  g_mst <- mst(g, algorithm = "prim")
  
  # Compute pairwise shortest paths over the MST
  sp <- distances(g_mst)
  
  # Project the shortest-path matrix to d dimensions
  res <- svds(sp, k = d)
  L <- diag(res$d, nrow = d, ncol = d)
  V <- res$v
  
  return(t(sqrt(L) %*% t(V)))
}

c_score <- function(ref){
  err <- 0
  for(i in 1:(nrow(ref)-1)){
    if(ref$diseaseType[i] != ref$diseaseType[i + 1]){
      err <- err + 1
    }
  }
  # Remove the actual switches of disease types from the error
  err <- err - (length(unique(ref$diseaseType)) - 1)
  c_scr <- (nrow(ref) - err)/nrow(ref)
  return(c_scr)
}

load("../data/DisGeNETv5.RData")
load("../results/separation.RData")

# Get the disease names and their types
dg <- disgenet %>% 
  select(diseaseName, diseaseType) %>% 
  filter(!duplicated(diseaseName)) %>% 
  arrange(diseaseName) %>% 
  mutate(diseaseType = as.factor(str_wrap(diseaseType, 50)))

# Colours for the disease types
cp <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
        "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#dfc27d", "#b15928", 
        "#878787")

# Rescale the disease separation matrices to have only positive values
D <- dH.matrix - min(dH.matrix)
diag(D) <- 0
Ds <- dS.matrix - min(dS.matrix)
diag(Ds) <- 0

# Perform non-centred Minimum Curvilinear Embedding (ncMCE)
mceH <- ncmce(D, d = 2)
mceS <- ncmce(Ds, d = 2)
mceJ <- ncmce(dJ.matrix, d = 2)

dg <- dg %>% mutate(mH1 = mceH[, 1], mH2 = mceH[, 2],
                    mS1 = mceS[, 1], mS2 = mceS[, 2],
                    mJ1 = mceJ[, 1], mJ2 = mceJ[, 2])

pH <- ggplot(dg, aes(mH1, mH2, colour = diseaseType)) + 
  geom_point() + scale_colour_manual(values = cp) + 
  labs(x = "Dim 1", y = "Dim 2") +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

pS <- ggplot(dg, aes(mS1, mS2, colour = diseaseType)) + 
  geom_point() + scale_colour_manual(values = cp) + 
  labs(x = "Dim 1", y = "Dim 2") +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

pJ <- ggplot(dg, aes(mJ1, mJ2, colour = diseaseType)) + 
  geom_point() + scale_colour_manual(values = cp) + 
  labs(x = "Dim 1", y = "Dim 2") +
  theme_bw() + theme(legend.position = "right", legend.title = element_blank())

# Evaluate the quality of the clustering on Dim 2
cH <- c_score(arrange(dg, mH2))
cS <- c_score(arrange(dg, mS2))
cJ <- c_score(arrange(dg, mJ2))

save(dg, pH, pS, pJ, cp, cH, cS, cJ, file = "../results/ndr_separation.RData")
