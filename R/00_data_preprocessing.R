options(stringsAsFactors = FALSE)
library("igraph")
library(org.Hs.eg.db)
library(biomaRt)

#############################################################
#
#       GWAS+OMIM gene-disease associations
#
#############################################################

#Read data and get rid of useless info.
menche <- read.table("data/DataS2_Menche_disease_genes.tsv", header = T, sep = "\t", quote = "", stringsAsFactors = F)
menche <- menche[,c("disease", "OMIM_genes", "GWAS_genes")]

#Combine OMIM and GWAS associations
#Split disease genes and put data in data frame
for (i in 1:nrow(menche)) {
  menche$disease.genes[i] <- lapply(i, function(x) unique(unlist(strsplit(paste(menche$OMIM_genes[i], menche$GWAS_genes[i], sep = ";"), ";"))))
}

menche <- menche[,c("disease", "disease.genes")]

#Save data frame
save(menche, file = "data/menche.RData")

#Now filter out redundant diseases (manual inspection)
#Row numbers to remove
to.remove <- c(1, 11, 12, 19, 260, 22, 31, 35, 41, 47, 50:52, 58, 63, 66, 70, 71, 74, 78, 79, 83, 85, 89, 97, 
               103, 104, 106, 108, 125, 126, 136, 141, 149, 153, 156:160, 163, 165, 166, 169, 170, 173, 174, 177, 
               184, 186, 193:195, 201, 205, 213, 230, 235, 237:239, 246, 255, 256, 258, 259, 266:269, 276, 
               278, 281, 290, 291, 293, 294, 295, 297, 299)

menche <- menche[-to.remove, ]

save(menche, file = "data/menche_filtered.RData")

#############################################################
#
#             HIPPIE network construction
#
#############################################################

#Read data and get rid of useless info.
hippie <- read.table("data/hippie_current.txt", header = FALSE, sep = "\t", quote = "", col.names = c("UniprotA", "EntrezA", "UniprotB", "EntrezB", "weight", "Desc"))
hippie <- hippie[, c("EntrezA", "EntrezB", "weight")]
# Search for empty entries and commas, didn't find any
sum(hippie$EntrezA == "")
length(grep(",", hippie$EntrezA))
sum(hippie$EntrezB == "")
length(grep(",", hippie$EntrezB)) 
# remove duplicates
hippie <- hippie[!duplicated(hippie), ]

#Apply filter to interactions (median confidence score)
qt <- quantile(hippie$weight)
# get weights above median
hippie.med <- hippie[hippie$weight >= qt["50%"],]

#Create igraph
gMedian <- graph_from_data_frame(hippie.med, directed = FALSE)
gMedian <- simplify(gMedian, edge.attr.comb = "min")

#Map Entrez gene IDs to Symbols for the median network (results are in same order as in gMedian)
map <- select(org.Hs.eg.db, keys = V(gMedian)$name, keytype = "ENTREZID", columns = c("ENTREZID", "SYMBOL"))
#Manually obtained symbols
map$SYMBOL[which(is.na(map$SYMBOL))] <- c("OCLN", "EPR1", "KCNJ12", "LOC100290070", "FSIP2", "TIMM23B", "ARHGAP27", "LOC100290337", "DUSP22", 
                                          "XAGE2", "CCDC7", "AREG", "PRRT2", "ITGAV", "COP", "YWHAB", "CASP12", "LOC653877", "GUSBP1", "RCC1", "XAGE1E", 
                                          "XAGE1B", "XAGE1E", "AE01", "RPS21", "LOC730426", "HSFX1", "ANXA8L1", "LOC100293069", "LOC100293130", "NANOG", 
                                          "LOC100509457", "MICA", "RPS17", "TRIM49C", "FRMPD2B", "FAM127B", "LOC100293351", "GOLGA6L9", "PPIAL4A", "PNMA6A", 
                                          "PNMA6B", "XK", "HIST2H3B", "UGP2", "LOC440264", "IGH", "LOC100293737", "DLGAP1", "SSX10", "UBE2L1", "SCX", 
                                          "LOC170549", "IGHV3OR15-7", "LOC90462", "LOC100293977", "TXLNGY", "MAGOH2P", "ZNF658B", "SPANXB1", "SPANXB1", 
                                          "WASHC2A", "HGH1", "LOC100290651", "HERC2P4", "LOC100292202", "ARL17B", "TSHZ2", "TSHZ2", "CFAP47", "AGAP10P", "TIMM23B")

V(gMedian)$symbol <- map$SYMBOL

# update Entrez genes IDs
V(gMedian)$name[V(gMedian)$name == 4950] <- 100506658
V(gMedian)$name[V(gMedian)$name == 727738] <- 374
V(gMedian)$name[V(gMedian)$name == 7449] <- 3685
V(gMedian)$name[V(gMedian)$name == 751867] <- 1104
V(gMedian)$name[V(gMedian)$name == 100293888] <- 79923
V(gMedian)$name[V(gMedian)$name == 3509] <- 3492
V(gMedian)$name[V(gMedian)$name == 50818] <- 112476
V(gMedian)$name[V(gMedian)$name == 4663] <- 7504

#Finally, extract the largest connected component
co <- clusters(gMedian)
gMedian <- induced_subgraph(gMedian, vids = which(co$membership == which.max(co$csize)))

#Save igraph and tab-separated files
save(gMedian, file = "data/gMedianNetwork.RData")
write.table(as_edgelist(gMedian), file = "data/gMedianNetwork.tsv", sep = "\t", quote = F, row.names = F, col.names = c("EntrezA", "EntrezB"))
write.table(as_edgelist(gMedian, names = F), file = "data/gMedianNetwork_forEmbeddding.tsv", sep = " ", quote = F, row.names = F, col.names = F)

write.table(as_edgelist(hPIN), file = "mogon_scripts/DIAMOnD-master/PPI.txt", sep = ",", quote = F, row.names = F, col.names = F)
