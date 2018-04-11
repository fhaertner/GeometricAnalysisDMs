options(stringsAsFactors = F)
library(cowplot)
library(igraph)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(tsne)
library(purrr)
load("data/hPIN_150k.RData")
load("data/coords_hPIN_150k.RData")
load("data/DisGeNETv5.RData")
load("data/hPIN_150k_H.RData")

# Given two sets of genes associated with a disease, computes the Jaccard index
# between the sets
jaccard <- function(S1, S2){
  return(length(intersect(S1, S2))/length(union(S1, S2)))
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
  return(mean.ij - mean(c(mean.ii, mean.jj)))
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

# retrieve the diseases and their associated proteins
diseases <- disgenet %>% 
  split(.$diseaseName)
# shorten the names to make a more readable plot later
names(diseases) <- shorten_terms(names(diseases))

# Pairwise separation based on hyperbolic distances
dH.matrix <- matrix(data = 0, nrow = length(diseases), ncol = length(diseases),
                      dimnames = list(names(diseases), names(diseases)))

# Pairwise separation based on shortest-paths
dS.matrix <- dH.matrix

# Compute pairwise shortest-paths
sp <- distances(hPIN, algorithm = "unweighted")

# Compute Jaccard similarity between disorders
by_dis <- disgenet %>% 
  split(.$diseaseName) %>% 
  map(~as.integer(.$geneId))
jc <- map(by_dis, function(x) map(by_dis, function(y) jaccard(x, y))) %>% 
  unlist() %>% 
  matrix(length(by_dis), length(by_dis), 
         byrow = TRUE, 
         dimnames = list(names(by_dis), names(by_dis)))
dJ.matrix <- 1 - jc

# compute the pairwise distances between the diseases and store in dH.matrix
for(i in 1:length(diseases)) {
  for(j in i:length(diseases)) {
    if(i != j) {
      idx.i <- which(V(hPIN)$name %in% diseases[[names(diseases)[i]]]$geneId)
      idx.j <- which(V(hPIN)$name %in% diseases[[names(diseases)[j]]]$geneId)
      d.ij <- get_distance(H, idx.i, idx.j)
      dH.matrix[i, j] <- d.ij
      dH.matrix[j, i] <- d.ij
      
      d.ij <- get_distance(sp, idx.i, idx.j)
      dS.matrix[i, j] <- d.ij
      dS.matrix[j, i] <- d.ij
    }
  }
}

# normalize the values in the matrix
dH.matrix.n <- transform_data(dH.matrix)
diag(dH.matrix.n) <- 0
dendr <- as.dendrogram(hclust(as.dist(dH.matrix.n), method = "ward.D"))

# get the diseases and the types they belong to
types <- data.frame(diseaseName = disgenet$diseaseName,
                    diseaseType = disgenet$diseaseType)
types <- types[!duplicated(types$diseaseName),]
types <- types[order(types$diseaseName),]
colorCodes <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
                "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#dfc27d", "#b15928", 
                "#878787")
names(colorCodes) <- sort(unique(types$diseaseType))


# extract a data frame from the dendrogram
off_raise <- 0
off_flat <- 70
ddata <- dendr %>% as.dendrogram %>% flatten.dendrogram(new_height = off_flat)
attr(ddata, "height") <- 75 
ddata <- ddata %>% raise.dendrogram(off_raise) %>% as.ggdend(type = "rectangle")
# ddata <- dendro_data(dendr, type = "rectangle")
# find the right order to assign the disease names their types and colour
types$order <- order(ddata$labels$label)
types <- types[order(types$order),]
groupCodes <- types$diseaseType

# this is for the p-value annotations
dend <- result %>% as.dendrogram %>% raise.dendrogram(off_raise) %>% flatten.dendrogram(new_height = off_flat) %>% hang.dendrogram
xy <- dend %>% get_nodes_xy
is_internal_node <- is.na(dend %>% get_nodes_attr("leaf"))
is_internal_node[which.max(xy[,2])] <- FALSE
xy <- xy[is_internal_node,]
xy <- xy[order(xy[,2]),]
nodes_range <- result$hclust$merge
node_pv <- result$edges$bp[-length(result$edges$bp)]

ddata$segments$y <- ddata$segments$y/5
ddata$segments$yend <- ddata$segments$yend/5
ddata$labels$label <- gsub("\n", "", ddata$labels$label)
# ddata$labels$label <- gsub("\\|.*","",ddata$labels$label)
p <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) +
  geom_text(data = ddata$labels, 
            aes(x = x, y = y, label = label,
                colour=types$diseaseType),
            size = 3.5, vjust = 0.5, hjust=0.001) +
  labs(color = element_blank()) + 
  annotate("text", x = xy[,1]+0.4, y = xy[,2]/5+0.2, label = round(node_pv, 2)*100, size = 4) +
  annotate("text", x = xy[155,1]+0.4, y = xy[155,2]/2+0.7, label = round(node_pv[155], 2)*100, size = 4) +
  theme(legend.position = "none") +
  scale_color_manual(values = colorCodes) +
  theme_dendro()
p

save_plot(filename = paste0("figs/ggdendro_plot.pdf"), plot = p, nrow = 1,
          ncol = 1, base_width = 37, base_height = 30, limitsize=F)
save(dH.matrix, dJ.matrix, dS.matrix, dendr, result, file = "results/separation.RData")
