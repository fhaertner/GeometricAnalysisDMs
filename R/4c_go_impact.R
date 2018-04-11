library(dplyr)
library(igraph)
library(FunEnrich)

load("data/hPIN_150k.RData")
load("data/DisGeNETv5.RData")
load("results/nav_impact_50epochs.RData")

imp <- mutate(imp, disease = names(eff))

# Get the list of proteins of DMs with significant impact
genes_imp <- unique(disgenet$geneId[disgenet$diseaseName %in% 
                               imp$disease[imp$pval_eff < 0.05]])

genes_nim <- unique(disgenet$geneId[disgenet$diseaseName %in% 
                               imp$disease[imp$pval_eff >= 0.05]])

# Check which of these proteins are in the hPIN
genes_imp <- genes_imp[genes_imp %in% V(hPIN)$name]
genes_nim <- genes_nim[genes_nim %in% V(hPIN)$name]

# Perform the enrichment analysis using the network as background
enr_imp <- fun_enrich(as.character(genes_imp), V(hPIN)$name, benjamini = TRUE)
enr_nim <- fun_enrich(as.character(genes_nim), V(hPIN)$name, benjamini = TRUE)

p_imp <- plot_fun_enrich(enr_imp, benjamini = TRUE)
p_nim <- plot_fun_enrich(enr_nim, benjamini = TRUE)

save(enr_imp, enr_nim, p_imp, p_nim, file = "results/go_impact.RData")