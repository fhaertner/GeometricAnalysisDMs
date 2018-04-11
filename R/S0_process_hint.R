library(igraph)
library(dplyr)
library(readr)

library(RCurl)
library(jsonlite)

hint.db <- read_tsv("../data/HINT_HomoSapiens_binary_hq_190218.txt", quote = "")
hint.db <- hint.db %>% 
  filter(!is.na(Alias_A) & !is.na(Alias_B)) %>% 
  select(Gene_A, Gene_B)

nodes <- union(hint.db$Gene_A, hint.db$Gene_B)
nodes <- sapply(strsplit(nodes, "[|]"), function(x) x[[1]][1])
nodes <- tibble(old = nodes, entrez = character(length(nodes)), 
                symbol = character(length(nodes)), 
                uniprot = character(length(nodes)))

coding <- character(length = nrow(nodes))
for(i in 1:nrow(nodes)){
  res <- tryCatch({
    gene.data <- fromJSON(getURL(paste0("http://mygene.info/v3/query?q=", 
                                        nodes$old[i],
                                        "&fields=entrezgene,symbol,uniprot,type_of_gene&species=human", 
                                        sep = "")))
    nodes$entrez[i] <- as.character(gene.data$hits$entrezgene)
    nodes$symbol[i] <- gene.data$hits$symbol
    nodes$uniprot[i] <- ifelse(is.null(gene.data$hits$uniprot$`Swiss-Prot`), "", 
                               gene.data$hits$uniprot$`Swiss-Prot`[1])
    coding[i] <- gene.data$hits$type_of_gene
  }, error = function(err){
    print(paste0("Symbol ", nodes$old[i], 
                 " most likely withdrawn from NCBI gene..."))
    return(NA)
  })
  if(is.na(res)){
    nodes$entrez[i] <- ""
  }
}
nodes$uniprot <- sapply(nodes$uniprot, function(x) x[1])
nodes[5014, 2:4] <- c("1589", "CYP21A2", "P08686")
nodes[6247, 2:4] <- c("64396", "GMCL1P1", "Q8NEA9")
node <- filter(nodes, entrez != "")
  
# Update node Entrez ID
hint.db$Gene_A <- left_join(hint.db, nodes, by = c("Gene_A" = "old"))$entrez
hint.db$Gene_B <- left_join(hint.db, nodes, by = c("Gene_B" = "old"))$entrez
hint.db <- filter(hint.db, Gene_A != "" & Gene_B != "")
nodes <- select(nodes, -old)
coding <- coding[!duplicated(nodes$entrez)]
nodes <- filter(nodes, !duplicated(nodes$entrez))

#Create an igraph object with the network data
hint <- graph_from_data_frame(hint.db, directed = F, vertices = nodes)

#Get rid of vertices with name = "" and non-coding genes and simplify
hint <- delete_vertices(hint, v = which(V(hint)$name == ""))
hint <- delete_vertices(hint, v = which(coding != "protein-coding"))
hint <- simplify(hint, remove.multiple = T, remove.loops = T)

# Consider the LCC only
cmp <- components(hint)
hint <- induced_subgraph(hint, vids = which(cmp$membership == 
                                              which.max(cmp$csize)))

saveRDS(hint, file = "../data/hint.rds")
