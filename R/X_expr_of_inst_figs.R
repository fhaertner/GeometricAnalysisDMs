
library(cowplot)
library(ggridges)
library(dplyr)
library(stringr)

# Part A ------------------------------------------------------------------

load("../results/dm_stats.RData")

partA <- ggplot(dm_stats, aes(lcc/dm_size, -log10(pval_h2))) + geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "red") +
  geom_vline(xintercept = 1, linetype = 2, colour = "blue") +
  labs(x = "Fraction of proteins in LCC", y = expression(-log[10](p-value))) +
  theme_bw()

# Part B ------------------------------------------------------------------

load("../data/coords_hPIN_150k.RData")
load("../data/DisGeNETv5.RData")

disgenet <- disgenet %>% 
  filter(diseaseName %in% c("cystic fibrosis", "fatty liver|steatohepatitis",
                            "myeloid leukemia", "polycystic ovary syndrome",
                            "small cell carcinoma of lung"))
coords$id <- as.integer(coords$id)
gene_coords <- left_join(disgenet, coords, by = c("geneId" = "id"))
gene_coords$diseaseName <- str_wrap(gene_coords$diseaseName, width = 15)

partB <- ggplot(gene_coords, aes(x = theta, y = diseaseName)) + 
  geom_density_ridges(bandwidth = 0.35) + 
  labs(x = expression(theta), y = "") + theme_bw()

# Part C ------------------------------------------------------------------

load("../results/nav_impact.RData")

partC <- ggplot(imp, aes(x = imp_eff_dis, y = imp_eff_rnd)) +
  geom_point(aes(colour = -log10(pval_eff), shape = (pval_eff < 0.05))) +
  coord_cartesian(xlim = c(-0.1, 0)) +
  labs(x = "Impact of disease proteins", y = "Impact of random proteins") + 
  scale_colour_gradient(low = "#fdb462", high = "blue", name = "-log(p-value)") +
  scale_shape(breaks = c(T, F), labels = c("significant", "not significant"), 
              name = element_blank()) + 
  #scale_size(name = "-log(p-value)") +
  theme_bw()

save_plot(filename = "figs/nav_impact.pdf", plot = partC,
          nrow = 1, ncol = 1, base_aspect_ratio = 2.2)

fig <- ggdraw() + 
  draw_plot(partA, x = 0, y = 0.5, width = 0.5, height = 0.5) + 
  draw_plot(partB, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot(partC, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(letters[1:3], x = c(0, 0.5, 0), y = c(1, 1, 0.5), size = 15)

save_plot("../figs/ei.pdf", fig, ncol = 2, nrow = 2)


# enrichment analysis
options(stringsAsFactors = F)
library(igraph)
library(reshape2)
load("results/nav_impact.RData")
load("data/protein_classes.RData")
load("data/hPIN_150k.RData")

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
  contingency <- matrix(0, nrow = 2, ncol = 2)
  contingency[1, 1] <- gene.counts["TRUE"] # k
  contingency[1, 2] <- bg.counts["TRUE"] # K - k
  contingency[2, 1] <- gene.counts["FALSE"] # n - k
  contingency[2, 2] <- bg.counts["FALSE"] # N + k - n - K
  # remove 'NA's if there are any
  contingency[which(is.na(contingency))] <- 0
  pval <- fisher.test(contingency, alternative = "greater")$p.value
  return(pval)
}

protein.names <- c("symbol", "tf", "rec", "trans")
protein.classes <- protein.classes[protein.names]
protein.names <- protein.names[-1]
# get background
bg <- lapply(eff, function(x) unlist(x$faulty_dis))
bg <- unname(unlist(bg))
# get the symbols from the id's through the hPIN
gene_ids <- unique(bg)
# create map
map_to_symbol <-V(hPIN)$symbol[which(V(hPIN)$name %in% gene_ids)] 
names(map_to_symbol) <- V(hPIN)$name[which(V(hPIN)$name %in% gene_ids)]
bg <- unname(map_to_symbol[as.character(bg)])
# count each gene in background
bg <- protein.classes[which(protein.classes$symbol %in% bg),]

disease.pvals <- data.frame(disease = names(eff), tf = rep(0, length(eff)), rec = rep(0, length(eff)), trans = rep(0, length(eff)), pval_eff = imp$pval_eff)
for (i in 1:length(eff)) {
  genes <- unlist(eff[[i]]$faulty_dis)
  genes <- unname(map_to_symbol[as.character(genes)])
  genes <- protein.classes[which(protein.classes$symbol %in% genes),]
  p.vals <- lapply(protein.names, function(x) do_fisher_tests(table(genes[x]), table(bg[x])))
  p.vals <- unlist(p.vals)
  names(p.vals) <- protein.names
  disease.pvals$tf[i] <- p.vals["tf"]
  disease.pvals$trans[i] <- p.vals["trans"]
  disease.pvals$rec[i] <- p.vals["rec"]
}



# do the same for the diseases
load("data/hPIN_150k.RData")
load("data/DisGeNETv5.RData")

diseases <- disgenet %>% 
  split(.$diseaseName)

# get background
bg <- lapply(diseases, function(x) unlist(x$geneSymbol))
bg <- unname(unlist(bg))
# count each gene in background
bg <- protein.classes[which(protein.classes$symbol %in% bg),]
disease.pvals <- data.frame(disease = names(diseases), tf = rep(0, length(diseases)), rec = rep(0, length(diseases)), trans = rep(0, length(diseases)), pval_eff = imp$pval_eff)
for (i in 1:length(diseases)) {
  genes <- unlist(diseases[[i]]$geneSymbol)
  genes <- protein.classes[which(protein.classes$symbol %in% genes),]
  p.vals <- lapply(protein.names, function(x) do_fisher_tests(table(genes[x]), table(bg[x])))
  p.vals <- unlist(p.vals)
  names(p.vals) <- protein.names
  disease.pvals$tf[i] <- p.vals["tf"]
  disease.pvals$trans[i] <- p.vals["trans"]
  disease.pvals$rec[i] <- p.vals["rec"]
}

# do plot
df.plot <- melt(disease.pvals, id = c("disease", "pval_eff"))
ea.plot <- ggplot(data = df.plot, aes(x=variable, y=value)) + geom_jitter(aes(colour = (pval_eff < 0.05))) +
  scale_colour_manual(breaks = NULL, labels = NULL, 
                      name = NULL, values = c("#fdb462", "blue")) +
  labs(x = "Protein class", y = "p-value") +
  theme_bw()


# count protein classes
gene_ids <- lapply(eff, function(x) unlist(x$faulty_dis))
gene_symbols <- lapply(gene_ids, function(x) unname(map_to_symbol[as.character(x)]))
classes.num <- data.frame(disease = names(gene_symbols), tf = rep(0, length(gene_symbols)), rec = rep(0, length(gene_symbols)), trans = rep(0, length(gene_symbols)), pval_eff = imp$pval_eff)
for (i in 1:length(gene_symbols)) {
  genes <- unlist(gene_symbols[i])
  genes <- protein.classes[which(protein.classes$symbol %in% genes),]
  classes.num$tf[i] <- sum(genes$tf == T)/nrow(genes)
  classes.num$trans[i] <- sum(genes$rec == T)/nrow(genes)
  classes.num$rec[i] <- sum(genes$trans == T)/nrow(genes)
}

# do plot
df.plot <- melt(classes.num, id = c("disease", "pval_eff"))
relative.plot <- ggplot(data = df.plot, aes(x=variable, y=value)) + geom_jitter(aes(colour = (pval_eff < 0.05))) +
  scale_colour_manual(breaks = c(T, F), labels = c("significant", "not significant"), 
                      name = "Impact of faulty proteins", values = c("#fdb462", "blue")) +
  labs(x = "Protein class", y = "Relative number of class proteins") +
  theme_bw()


final.plot <- plot_grid(ea.plot, relative.plot, labels = c("A", "B"), rel_widths = c(2, 3))
save_plot("figs/protein_classes_faulty.pdf", final.plot, ncol = 2, nrow = 1, base_aspect_ratio = 1.3)
