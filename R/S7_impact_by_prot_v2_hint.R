library(dplyr)
library(igraph)
library(cowplot)
library(FunEnrich)

load("results/S_nav_hint_50epochs.RData")
load("data/protein_classes.RData")
hint_net <- readRDS("data/hint.rds")

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
freq_faul <- prot_counts[which(prot_counts > quantile(prot_counts)["75%"])]
infreq_faul <- prot_counts[which(prot_counts <= quantile(prot_counts)["75%"])]

enr_freq <- fun_enrich(names(freq_faul), V(hint_net)$name, benjamini = TRUE)
p_freq <- plot_fun_enrich(enr_freq, benjamini = TRUE, char_per_line = 40)

enr_infreq <- fun_enrich(names(infreq_faul), V(hint_net)$name, benjamini = TRUE)
p_infreq <- plot_fun_enrich(enr_infreq, benjamini = TRUE, char_per_line = 40)

sym_freq <- V(hint_net)$symbol[V(hint_net)$name %in% names(freq_faul)]
sym_infreq <- V(hint_net)$symbol[V(hint_net)$name %in% names(infreq_faul)]

df <- tibble(status = rep(c("FDA-approved", "Potential drug target"), 2),
             faulty = rep(c("Frequent faulty", "Infrequent faulty"), each = 2),
             val = c(sum(sym_freq %in% protein.classes$symbol[protein.classes$fda]),
                      sum(sym_freq %in% protein.classes$symbol[protein.classes$pot_targets]),
                      sum(sym_infreq %in% protein.classes$symbol[protein.classes$fda]),
                      sum(sym_infreq %in% protein.classes$symbol[protein.classes$pot_targets])),
             tot = c(rep(length(sym_freq), 2), rep(length(sym_infreq), 2)))

p_dt <- ggplot(df, aes(faulty, val/tot, fill = status)) + 
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#ef8a62", "#67a9cf")) +
  labs(x = "", y = "Fraction", fill = "") + theme_bw() + 
  theme(legend.background = element_blank(), 
        legend.justification = c(1,1), legend.position = c(1,1))

# Fisher's test to ensure the proportions are significant
m_fda <- matrix(c(sum(sym_freq %in% protein.classes$symbol[protein.classes$fda]),
                   sum(sym_infreq %in% protein.classes$symbol[protein.classes$fda]),
                   sum(!(sym_freq %in% protein.classes$symbol[protein.classes$fda])),
                   sum(!(sym_infreq %in% protein.classes$symbol[protein.classes$fda]))), 
                 2, 2, byrow = TRUE)
pval_fda <- fisher.test(m_fda, alternative = "greater")$p.value


m_pot <- matrix(c(sum(sym_freq %in% protein.classes$symbol[protein.classes$pot_targets]),
                  sum(sym_infreq %in% protein.classes$symbol[protein.classes$pot_targets]),
                  sum(!(sym_freq %in% protein.classes$symbol[protein.classes$pot_targets])),
                  sum(!(sym_infreq %in% protein.classes$symbol[protein.classes$pot_targets]))), 
                2, 2, byrow = TRUE)
pval_pot <- fisher.test(m_pot, alternative = "greater")$p.value

freq_tab <- tibble(entrez = names(freq_faul), times_used_as_faulty = freq_faul, 
                   fda = sym_freq %in% protein.classes$symbol[protein.classes$fda],
                   potential_drug_target = sym_freq %in% protein.classes$symbol[protein.classes$pot_targets])
freq_tab <- arrange(freq_tab, desc(times_used_as_faulty))

save(df, p_dt, pval_fda, pval_pot, freq_tab, 
     enr_freq, p_freq, enr_infreq, p_infreq, file = "results/S_nav_imp_drugs_hint.RData")
