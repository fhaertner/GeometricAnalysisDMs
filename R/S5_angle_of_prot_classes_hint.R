library(cowplot)
library(igraph)
library(ggridges)
library(dplyr)
library(stringr)

#' Given a network and its hyperbolic coordinates, get the angular components
#' of its proteins and the classes they belong to
#' 
#'  @param net igraph object object representing the hPIN
#'  @param coords data frame with columns 'id' for protein id and 'theta' for
#'  the angular component
#'  @param classes data frame that maps protein symbols to a class
#'
get_angle_distr <- function(net, coords, classes) {
  class.names <- colnames(classes[, -1])
  
  thetas <- c()
  types <- c()
  for(i in 1:length(class.names)){
    # Determine which network proteins are in class
    idx <- which(V(net)$symbol %in% classes$symbol[classes[, class.names[i]]])
    thetas <- c(thetas, coords$theta[idx])
    types <- c(types, rep(class.names[i], length(idx)))
  }
  df <- data.frame(theta = thetas, class = factor(types,
                                                  levels = unique(types),
                                                  ordered = T))
  return(df)
}

#' Given a network and its hyperbolic coordinates, plot the angular
#' distribution of its protein classes
#' 
#'  @param net igraph object object representing the hPIN
#'  @param coords data frame with columns 'id' for protein id and 'theta' for
#'  the angular component
#'  @param classes data frame that maps protein symbols to a class
#'  @param bins integer for the number of bins
#'
plot_theta_distr <- function(net, coords, classes, bins){
  df <- get_angle_distr(net, coords, classes)
  p <- ggplot(df, aes(theta, ..density.., fill = class)) + 
    geom_histogram(position = "identity", bins = bins) + 
    scale_fill_manual(values = c("#ED5051", "#FE880F", "#45A939", "#FB9A99",
                                 "#87C196", "#A6CEE3")) +
    labs(x = expression(theta), y = "Density") + theme_bw() + 
    theme(legend.title = element_blank(), legend.background = element_blank(),
          legend.justification = c(0,1), legend.position = c(0,1))
  return(p)
}

#' Perform Fisher's test
#' 
#' Carry out the required over-representation tests for each given term
#' 
#' @param gene.counts The term counts for the genes of interest
#' @param bg.counts The term counts for the background
#' 
#' @return A vector with the p-values for each term
#' 
do_fisher_tests <- function(gene.counts, bg.counts){
  
  pvals <- numeric(length = length(gene.counts))
  
  contingency <- matrix(0, nrow = 2, ncol = 2)
  j <- 1
  for(i in names(gene.counts)){
    contingency[1, 1] <- gene.counts[i] # k
    contingency[1, 2] <- bg.counts[i] - gene.counts[i] # K - k
    contingency[2, 1] <- sum(gene.counts) - gene.counts[i] # n - k
    contingency[2, 2] <- sum(bg.counts) + gene.counts[i] - sum(gene.counts) -
      bg.counts[i] # N + k - n - K
    pvals[j] <- fisher.test(contingency, alternative = "greater")$p.value
    j <- j+1
  }
  return(pvals)
}

#' Bin the proteins and count frequency of each class. Has two modes:
#' 'frequency' just determines the most frequent class for each bin
#' 'enrichment' does an enrichment analyses and determines the most
#' enriched class for each bin
#' 
#' 
#' @param protein_classes data frame storing which classes each protein
#'  belongs to
#' @param pin igraph object object representing the hPIN
#' @param coords data frame with columns 'id' for protein id and 'theta' for
#'  the angular component
#' @param num_bins integer for the number of bins
#' @param mode character, can be 'frequency' or 'enrichment' for the two
#'  modes
#' 
#' @return A vector with the p-values for each term
#' 
get_class_intervals <- function(protein_classes, pin, coords, num_bins = 10,
                                 mode = "frequency") {
  # get the distribution and bin the angles
  theta_distr <- get_angle_distr(pin, coords, protein_classes)
  bins <- seq(min(theta_distr$theta), max(theta_distr$theta),
              (max(theta_distr$theta) - min(theta_distr$theta))/(num_bins - 1))
  theta_distr$bins <- findInterval(theta_distr$theta, bins)
  # create data frame storing most frequent class in each interval
  class_intervals <- data.frame(interval = 1:num_bins,
                                 class = rep(0, num_bins), range = bins)
  protein_names <- names(protein_classes[, -1])
  if(mode == "enrichment") {
    background <- protein_classes[which(protein_classes$symbol %in%
                                          V(pin)$symbol),]
    # go through the bins and count proteins of each class
    for (i in 1:nrow(class_intervals)) {
      window.counts <- table(theta_distr$class[
        theta_distr$bins == class_intervals$interval[i]])
      
      bg.counts <- c()
      for(j in 1:length(protein_names)){
        # Determine which network proteins are in class
        bg.counts <- c(bg.counts, rep(protein_names[j],
                                      sum(background[protein_names[j]] == T)))
      }
      bg.counts <- table(bg.counts)
      
      p.vals <- do_fisher_tests(window.counts, bg.counts)
      names(p.vals) <- names(window.counts)
      class_intervals$class[i] <-  names(which.min(p.vals))
    }
  }
  else if(mode == "frequency"){
    # go through the bins and determine which class occurs most often in bin
    for (i in 1:nrow(class_intervals)) {
      class_intervals$class[i] <-  protein_names[which.max(unlist(
        lapply(protein_names, function(n) length(grep(n, theta_distr$class[
          theta_distr$bins == class_intervals$interval[i]])))))]
    }
  }
  return(class_intervals) 
}

#' Generates a ggplot object displaying the protein distributions
#' 
#' 
#' @param class_intervals data frame storing the most frequent class for each
#'  bin
#' @param protein_names character list for the classes that could occur in the
#'  bins
#' @param coords character list for the diseases to be plotted
#' 
#' @return A vector with the p-values for each term
#' 
plot_ridges <- function(class_intervals, protein_names, diseases) {
  load("data/DisGeNETv5.RData")
  # remove rows with same consecutive classes
  to.remove <- c()
  for (i in 1:(nrow(class_intervals) - 1)) {
    if(class_intervals$class[i] == class_intervals$class[i + 1]) {
      to.remove <- c(to.remove, i)
    }
  }
  class_intervals <- class_intervals[-to.remove,]
  
  # select some diseases
  disgenet <- disgenet %>% 
    filter(diseaseName %in% diseases)
  coords$id <- as.integer(coords$id)
  gene_coords <- left_join(disgenet, coords, by = c("geneId" = "id"))
  gene_coords$diseaseName <- str_wrap(gene_coords$diseaseName, width = 15)
  gene_coords <- mutate(gene_coords, 
                        diseaseName = factor(diseaseName, 
                                             levels = str_wrap(diseases, width = 15), 
                                             ordered = TRUE))
  gene_coords <- filter(gene_coords, !is.na(theta))
  
  # for the background of the plot choose a color for each class
  colorCodes <- c("#ED5051", "#FE880F", "#45A939", "#FB9A99", 
                  "#87C196", "#A6CEE3")
  names(colorCodes) <- protein_names
  
  # generate plot
  final.plot <- ggplot(gene_coords, aes(x = theta, y = diseaseName)) + 
    geom_density_ridges(bandwidth = 0.1, alpha = 1) + 
    labs(x = expression(theta), y = "") + 
    annotate("rect", xmin = c(-Inf, class_intervals$range[-1]),
             xmax = c(class_intervals$range[-1], Inf),
             ymin = -Inf, ymax = Inf,
             alpha = 0.5, fill = colorCodes[class_intervals$class]) +
    annotate("text",
             x=(class_intervals$range + c(class_intervals$range[-1],
                                           2.2*pi))/2,
             y = rep(0.8, length(names(colorCodes[class_intervals$class]))),
             label = names(colorCodes[class_intervals$class])) +
    geom_density_ridges(bandwidth = 0.1, alpha = 1) + 
    theme_bw()
  return(final.plot)
}

# ==============================================
#
# Use bins to find most frequent protein classes
#
# ==============================================


# load the data
load("data/protein_classes.RData")
load("data/hint_coords.RData")
hint_net <- readRDS("data/hint.rds")
coords <- coords_lh$polar
# rename the columns in protein_classes
names(protein.classes) <- c("symbol", "TF", "Rec", "Trans", "Enzyme", "RBP",
                            "Skel", "Ubi", "Plasma", "Cancer", "FDA", "Pot_Targets")
protein.classes <- protein.classes[,-c(5, 9, 10, 11, 12)]
# class_intervals <- get_class_intervals(protein.classes, hPIN, coords,
#                                        mode = "frequency", num_bins = 22)
class_intervals <- get_class_intervals(protein.classes, hint_net, coords,
                                       mode = "enrichment", num_bins = 15)
p_supp <- plot_ridges(class_intervals, names(protein.classes[, -1]),
                           sort(c("myeloid leukemia", 
                             "myocardial infarction", "schizophrenia",
                             "celiac disease", "diabetes",
                             "alzheimer's disease", "rheumatoid arthritis"),
                             decreasing = TRUE))

save(p_supp, file = "results/dm_angles_hint.RData")

