options(stringsAsFactors = F)
library(pvclust)
library(igraph)
library(parallel)
load("data/hPIN_150k.RData")
load("data/DisGeNETv5.RData")
load("data/hPIN_150k_H.RData")

transform_data <- function(data) {
  min.point <- min(data)
  if(min.point < 0) {
    return(data - min.point)
  }
  return(data)
}

# shorten the labels if they are longer than 'cutoff' characters
shorten_terms <- function(terms, cutoff=50) {
  new_terms <- c()
  for(t in terms) {
    if(nchar(t) > 50) {
      # if the number of characters in the term is larger than 50, split it in
      # several lines, store the new lines in new_t
      new_t <- ""
      for(s in unlist(strsplit(t, " "))) {
        # check if current line has less than 50 characters if you add s, 
        # if yes, add s, otherwise start new line
        if((nchar(gsub(".*\n", "", new_t)) + nchar(s) + 1) < cutoff) {
          new_t <- paste(new_t, " ", s, sep = "")
        } else {
          new_t <- paste(new_t, "\n", s, sep = "")
        }
      }
      new_terms <- c(new_terms, trimws(new_t))
    } else {
      new_terms <- c(new_terms, t)
    }
  }
  return(new_terms)
}

#' Computes the distance between two sets of proteins in the hyperbolic plane
#' 
#'  @param H matrix, containing the pairwise hyperbolic distances between 
#'  proteins
#'  @param idx.1 indices in H of the the first set of proteins 
#'  @param idx.2 indices in H of the the second set of proteins 
#' 
get_distance <- function(H, idx.i, idx.j) {
  mean.ii <- get_mean_row_min_distance(H, idx.i, idx.i)
  mean.jj <- get_mean_row_min_distance(H, idx.j, idx.j)
  mean.ij <- get_mean_row_min_distance(H, idx.i, idx.j)
  result <- mean.ij - mean(c(mean.ii, mean.jj))
  return(result)
}

#' Computes the average minimum distance of the rows in a subsetted 
#' distance matrix
#' 
#'  @param H matrix, containing the pairwise hyperbolic distances between 
#'  proteins
#'  @param idx.1 indices in H of the the first set of proteins 
#'  @param idx.2 indices in H of the the second set of proteins 
#' 
get_mean_row_min_distance <- function(H, idx.i, idx.j) {
  # subset matrix
  sub.H <- H[idx.i, idx.j]
  diag(sub.H) <- Inf
  # get the minima of each row and compute the mean
  mean.ij <- mean(apply(sub.H, 1, FUN = function(x) min(x)))
  # subset matrix
  sub.H <- H[idx.j, idx.i]
  diag(sub.H) <- Inf
  # get the minima of each row and compute the mean
  mean.ji <- mean(apply(sub.H, 1, FUN = function(x) min(x)))
  return(mean(c(mean.ij, mean.ji)))
}

# retrieve the diseases and their associated proteins
diseases <- disgenet %>% 
  split(.$diseaseName)
max.genes <- unique(disgenet$geneId)

# construct data frame with coords of the associated disease proteins
dis.genes <- data.frame(matrix(
  rep(0, length(names(diseases))*length(max.genes)), 
  ncol=length(max.genes), dimnames = list(names(diseases), max.genes)))
names(dis.genes) <- max.genes
dis.genes <- cbind(dis.genes, data.frame(disease = names(diseases)))
for (i in 1:length(max.genes)) {
  disease.name <- disgenet$diseaseName[grep(paste("^", max.genes[i], 
                                                  "$", sep=""), 
                                            disgenet$geneId)]
  idx <- which(dis.genes$disease %in% disease.name)
  dis.genes[,i][idx] <- T
}
dis.genes <- t(dis.genes)
dis.genes <- as.data.frame(dis.genes)
dis.genes <- dis.genes[-length(dis.genes),]

compute_dist <- function(diseases) {
  dist.matrix <- matrix(data = rep(0, times = ncol(diseases)*ncol(diseases)),
                        ncol = ncol(diseases),
                        dimnames = list(diseases$disease, diseases$disease))
  gene_ids <- rownames(diseases)
  # compute the pairwise distances between the diseases and store in matrix
  for(i in 1:ncol(diseases)) {
    for(j in i:ncol(diseases)) {
      if(i != j) {
        disease_genes <- gene_ids[which(diseases[,i] == 1)]
        idx.i <- which(V(hPIN)$name %in% disease_genes)
        disease_genes <- gene_ids[which(diseases[,j] == 1)]
        idx.j <- which(V(hPIN)$name %in% disease_genes)
        d.ij <- get_distance(H, idx.i, idx.j)
        dist.matrix[i, j] <- d.ij
        dist.matrix[j, i] <- d.ij
      }
    }
  }
  # normalize the values in the matrix
  dist.matrix.n <- transform_data(dist.matrix)
  diag(dist.matrix.n) <- 0
  result <- as.dist(dist.matrix.n)
  attr(result, "method") <- "compute_dist"
  return(result)
}

# Calculate the number of cores
no_cores <- 8

# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, list("dis.genes", "H", "hPIN", "compute_dist",
                       "get_distance", "get_mean_row_min_distance", 
                       "transform_data", "diseases", "disgenet"))
clusterEvalQ(cl, list(library(igraph)))

result <- parPvclust(data=dis.genes, method.hclust = "ward.D", cl = cl, 
                     method.dist = compute_dist, nboot = 100)
save(result, file = "results/bootstrap_dendrogram.RData")