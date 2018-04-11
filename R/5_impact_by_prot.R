library(dplyr)
library(igraph)
library(reshape2)
library(cowplot)
library(FunEnrich)

options(stringsAsFactors = F)
load("results/nav_impact_50epochs.RData")
load("data/protein_classes.RData")
load("data/hPIN_150k.RData")
n <- length(eff[[1]]$faulty_dis[[1]])
epochs <- length(eff[[1]]$eff_dis)
prot_efficiencies <- data.frame()
for(i in 1:length(eff)) {
  eff_dis <- unlist(eff[[i]]$eff_dis)
  eff_dis <- as.vector(sapply(eff_dis, function(x) rep(x, n)))
  dis_prot <- unlist(eff[[i]]$faulty_dis)
  prot_efficiencies <- rbind(prot_efficiencies, data.frame(protein = dis_prot,
                                                           eff = eff_dis))
}

# find most frequently used proteins (e.g. proteins used more than 50 times)
prot_counts <- table(prot_efficiencies$protein)
max_prot <- prot_counts[which(prot_counts > quantile(prot_counts)["75%"])]

enr <- fun_enrich(names(max_prot), V(hPIN)$name, benjamini = TRUE)
p_max_prot <- plot_fun_enrich(enr, benjamini = TRUE, char_per_line = 40)




p_vals <- data.frame(protein = names(max_prot),
                     p_val = rep(0, length(max_prot)))
for(i in 1:nrow(p_vals)) {
  prot_eff <- prot_efficiencies %>%
    filter(protein == p_vals$protein[i]) %>%
    select(eff) %>%
    unlist()
  p_vals$p_val[i] <- wilcox.test(prot_eff, prot_efficiencies$eff,
                                 alternative = "less")$p.value
} 

# add information about proteins being fda drug targets or not
symbols_map <- data.frame(entrez = V(hPIN)$name, symbol = V(hPIN)$symbol)
p_vals <- left_join(p_vals, symbols_map, by = c("protein" = "entrez"))
protein.classes <- protein.classes %>%
  select(symbol, fda, rec, enzyme, plasma, pot_targets)
p_vals <- left_join(p_vals, protein.classes)

# see how many of the heavy impact proteins are also fda drug targets
sign_group <- p_vals %>% 
  filter(p_val < 0.05)
table(sign_group$fda)

# do enrichment analysis

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
  contingency[1, 2] <- bg.counts["TRUE"] - gene.counts["TRUE"] # K - k
  contingency[2, 1] <- sum(gene.counts) - gene.counts["TRUE"] # n - k
  contingency[2, 2] <- sum(bg.counts) + gene.counts["TRUE"] -
    sum(gene.counts) - bg.counts["TRUE"] # N + k - n - K
  pval <- fisher.test(contingency, alternative = "greater")$p.value
  return(pval)
}

bg_prot <- data.frame(protein = names(prot_counts))
bg_prot <- left_join(bg_prot, symbols_map, by = c("protein" = "entrez"))
bg_prot <- left_join(bg_prot, protein.classes)

classes <- c("fda", "rec", "enzyme", "plasma")
enrichment <- data.frame(class = classes, p_val = rep(0, length(classes)),
                         frac_sign = rep(0, length(classes)),
                         frac_bg = rep(0, length(classes)))
for(i in 1:nrow(enrichment)) {
  bg.counts <- table(bg_prot[enrichment$class[i]])
  gene.counts <- table(sign_group[enrichment$class[i]])
  enrichment$p_val[i] <- do_fisher_tests(gene.counts, bg.counts)
  enrichment$frac_sign[i] <- gene.counts["TRUE"]/sum(gene.counts)
  enrichment$frac_bg[i] <- bg.counts["TRUE"]/sum(bg.counts)
}
enrichment$class <- c("FDA drug target", "Rec", "Enzyme", "Plasma")

# do plot
df <- melt(enrichment,id.vars=c('class', 'p_val'), 
           measure.vars=factor(c('frac_sign','frac_bg')))

plot.impact <- ggplot(df, aes(class, value)) + 
  geom_bar(aes(fill = variable, colour=(p_val < 0.05)), 
           position = "dodge", stat="identity") +
  labs(y = "Fraction of class proteins", x = "Class") +
  scale_fill_discrete(limit = c('frac_sign','frac_bg'),
                   labels = c('Heavy impact proteins','All proteins'),
                   guide = guide_legend(title = NULL)) +
  scale_color_manual(values=c("black", "red"),
                     labels = c("Not enriched", "Enriched"),
                     guide = guide_legend(title = NULL)) + 
  theme_bw()

save_plot(plot.impact, base_aspect_ratio = 1.7, 
          filename = "figs/classes_vs_impact_50epochs.pdf")


sign_group <- sign_group[which(!is.na(sign_group$pot_targets)),]
pot_targets <- sign_group$protein[sign_group$pot_targets == T]

save(p_vals, sign_group, pot_targets, plot.impact, file = "impact_by_prot.RData")
