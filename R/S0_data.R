library(igraph)
library(dplyr)
library(stringr)


load("data/hPIN_150k.RData")
load("data/coords_hPIN_150k.RData")
load("data/protein_classes.RData")
# write hPIN to tsv file
write.table(as_edgelist(hPIN), file = "data/S_hPIN.tsv", sep = "\t", quote = F,
            row.names = F, col.names = c("EntrezA", "EntrezB"))

# gather information about proteins
protein_info <- data.frame(Entrez = V(hPIN)$name, Symbol = V(hPIN)$symbol,
                           Uniprot = V(hPIN)$uniprot)
protein_info <- left_join(protein_info, coords, by = c("Entrez" = "id"))
protein_info <- left_join(protein_info, protein.classes,
                          by = c("Symbol" = "symbol"))

# write to file
write.table(protein_info, file = "data/S_protein_info_hPIN.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)
rm(list = ls())


# disgenet
load("data/DisGeNETv5.RData")
write.table(disgenet, file = "data/S_DisGeNETv5.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)
rm(list = ls())


# hint network
load("data/hint_coords.RData")
load("data/protein_classes.RData")
hint_net <- readRDS("data/hint.rds")

# write hPIN to tsv file
write.table(as_edgelist(hint_net), file = "data/S_hint_pin.tsv", sep = "\t", quote = F,
            row.names = F, col.names = c("EntrezA", "EntrezB"))

# gather information about proteins
protein_info <- data.frame(Entrez = V(hint_net)$name, Symbol = V(hint_net)$symbol,
                           Uniprot = V(hint_net)$uniprot)
protein_info <- left_join(protein_info, coords_lh$polar, by = c("Entrez" = "id"))
protein_info <- left_join(protein_info, protein.classes,
                          by = c("Symbol" = "symbol"))

# write to file
write.table(protein_info, file = "data/S_protein_info_hint.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)
rm(list = ls())
