library(readr)
library(RCurl)
library(jsonlite)
library(dplyr)
library(NetHypGeom)

# Get current human PPI's from STRING database
string_main <- "https://stringdb-static.org/download/protein.links.v10.5/"
string_dataset <- "9606.protein.links.v10.5.txt.gz"
string_db <- read_table2(file = paste0(string_main, string_dataset))


# Get the mappings from string identifiers to entrez id's
map_main <- "https://string-db.org/mapping_files/entrez_mappings/"
map_dataset <- "entrez_gene_id.vs.string.v10.28042015.tsv"
map <- read_tsv(file = paste0(map_main, map_dataset),
                col_names = c("entrez","string_id"), skip = 1)

# mapping
string_db$entrezA <- left_join(string_db, map, 
                               by = c("protein1" = 
                                        "string_id"))$entrez
string_db$entrezB <- left_join(string_db, map, 
                               by = c("protein2" = 
                                        "string_id"))$entrez

nodes <- union(string_db$entrezA, string_db$entrezB)
nodes <- tibble(old = nodes, entrez = character(length(nodes)), 
                symbol = character(length(nodes)), 
                uniprot = character(length(nodes)))

#Some STRING entrez gene IDs might be out of date
#Lets update them using the mygene.info service
#In the process, let's get SYMBOL and UNIPROT information
coding <- character(length = nrow(nodes))
for(i in 1:nrow(nodes)){
  res <- tryCatch({
    gene.data <- fromJSON(
      getURL(paste0("http://mygene.info/v2/gene/", nodes$old[i],
                    "?fields=entrezgene,symbol,uniprot,type_of_gene",
                    sep = "")))
    nodes$entrez[i] <- as.character(gene.data$entrezgene)
    nodes$symbol[i] <- gene.data$symbol
    nodes$uniprot[i] <- ifelse(is.null(gene.data$uniprot$`Swiss-Prot`), "", 
                               gene.data$uniprot$`Swiss-Prot`)
    coding[i] <- gene.data$type_of_gene
  }, error = function(err){
    print(paste0("Entrez ", nodes$old[i], 
                 " most likely withdrawn from NCBI gene..."))
    return(NA)
  })
  if(is.na(res)){
    nodes$entrez[i] <- ""
  }
}

# Update node Entrez ID
string_db$entrezA <- left_join(string_db, nodes,
                               by = c("entrezA" = "old"))$entrez
string_db$entrezB <- left_join(string_db, nodes,
                               by = c("entrezB" = "old"))$entrez
nodes <- select(nodes, -old)
coding <- coding[!duplicated(nodes$entrez)]
nodes <- filter(nodes, !duplicated(nodes$entrez))

# remove entries whose mapping was not sucessfull
string_db <- string_db[-grep("^$", string_db$entrezA),]
string_db <- string_db[-grep("^$", string_db$entrezB),]

#Apply filter to get the highest scoring interactions
qt <- quantile(string_db$combined_score, c(.91, .92, .93))
# get weights above 92th percentile
string_db.filter <- string_db[string_db$combined_score >= qt["92%"],]
print(length(unique(string_db.filter$protein1)))

#Create igraph
string_pin <- graph_from_data_frame(string_db.filter, directed = FALSE)
string_pin <- simplify(string_pin, edge.attr.comb = "min")

#Finally, extract the largest connected component
co <- clusters(string_pin)
string_pin <- induced_subgraph(
  string_pin,
  vids = which(co$membership == which.max(co$csize)))

# determine the network temperature T
# begin by estimating the clustering at T=0
N <- vcount(string_pin)
avg.k <- mean(degree(string_pin))
gma <- fit_power_law(degree(string_pin))$alpha

for(i in 1:10){
  net <- ps_model(N = N, avg.k = avg.k, gma = gma, Temp = 0)
  clust.at.zero[i] <- transitivity(net$network, "average")
}
source("utility_functions.R")

# now estimate the network temperature

Temp <- estimate_temp(mean(clust.at.zero), transitivity(string_pin, "average"))

#Save igraph
save(string_pin, Temp, file = "data/StringNetwork.RData")
