
library(readr)
library(dplyr)
library(purrr)
library(tidyr)

# Given two sets of genes associated with a disease, computes the Jaccard index
# between the sets
jaccard <- function(S1, S2){
  return(length(intersect(S1, S2))/length(union(S1, S2)))
}

# Get gene-disease associations from DisGeNET v5.0
disgenet_main <- "http://www.disgenet.org/ds/DisGeNET/results/"
disgenet_dataset <- "all_gene_disease_pmid_associations.tsv.gz"
disgenet <- read_tsv(file = paste0(disgenet_main, disgenet_dataset), quote = "")

# Focus on diseaseType == 'disease' and associations supported by at least 3
# publications
disgenet <- disgenet %>% 
  filter(diseaseType == "disease") %>% 
  select(geneId, geneSymbol, diseaseId, diseaseName, pmid) %>%
  group_by(geneId, geneSymbol, diseaseId, diseaseName) %>% 
  summarise(NofPmids = n()) %>% 
  filter(NofPmids >= 3)

# Keep only diseases with >= 50 associated genes
disgenet <- disgenet %>% 
  group_by(diseaseId) %>% 
  mutate(assocGenes = n()) %>% 
  ungroup() %>% 
  filter(assocGenes >= 50) %>% 
  select(geneId, geneSymbol, diseaseName)

# Compute Jaccard index between diseases to detect synonyms (redundant info)
disgenet$diseaseName <- tolower(disgenet$diseaseName)

# Split data by disease
by_dis <- disgenet %>% 
  split(.$diseaseName) %>% 
  map(~as.integer(.$geneId))

# Compute Jaccard
jc <- map(by_dis, function(x) map(by_dis, function(y) jaccard(x, y))) %>% 
  unlist() %>% 
  matrix(length(by_dis), length(by_dis), 
         byrow = TRUE, 
         dimnames = list(names(by_dis), names(by_dis)))

# Cluster disorders with high gene-set overlap
clu <- hclust(as.dist(1 - jc))
# plot(clu, hang = -1) # To look at the dendrogram
clu_cut <- cutree(clu, h = 0.6)

# Merge redundant diseases by translating disease names into merged names
merge_dis <- function(dis_name){
  clu_id <- clu_cut[dis_name]
  to_merge <- names(clu_cut)[clu_cut == clu_id]
  return(paste(to_merge, collapse = "|"))
}

disgenet$diseaseName <- disgenet$diseaseName %>% map_chr(merge_dis)

# Finally, get rid of duplicates and bacterial and viral diseases
disgenet <- disgenet[!duplicated(disgenet), ]
bac_vir <- paste0("^sepsis|^helicobacter|^hepatitis|^human immuno|^influenza",
                  "|^malaria|^persistent embry|^pneum|^tubercu|^autosomal",
                  "|^graft")
idx <- grep(bac_vir, 
           disgenet$diseaseName)
disgenet <- disgenet[-idx, ]
# write(sort(unique(disgenet$diseaseName)), file = "../data/disease_types.tsv")

# Add information about disease types from revision 10 of the International 
# Statistical Classification of Diseases and Related Health Problems
types <- read_tsv("../data/disease_types.tsv", quote = "")
disgenet <- left_join(disgenet, types, by = "diseaseName")

save(disgenet, file = "../data/DisGeNETv5.RData")


