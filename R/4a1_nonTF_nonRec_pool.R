
library(igraph)
library(ggplot2)
library(dplyr)

load("../data/hPIN_150k.RData")
load("../data/protein_classes.RData")

#From the set of all proteins, form a pool that are not TF or Rec but that have 
# similar degrees
tf <- V(hPIN)$symbol %in% protein.classes$symbol[protein.classes$tf]
rec <- V(hPIN)$symbol %in% protein.classes$symbol[protein.classes$rec]

tf.deg <- degree(hPIN)[tf]
rec.deg <- degree(hPIN)[rec]

avail <- which(!(tf | rec))
not.tf <- vector("numeric", length = sum(tf))
not.rec <- vector("numeric", length = sum(rec))

for(i in 1:length(not.tf)){
  flag <- F
  curr.deg <- tf.deg[i]
  
  while(!flag){
    candidates <- degree(hPIN)[avail] == curr.deg
    if(sum(candidates) > 0){
      flag <- T
      not.tf[i] <- ifelse(sum(candidates) > 1, 
                          sample(avail[candidates], 1), 
                          avail[candidates])
    }else{
      curr.deg <- curr.deg - 1
    }
  }
}

for(i in 1:length(not.rec)){
  flag <- F
  curr.deg <- rec.deg[i]
  
  while(!flag){
    candidates <- degree(hPIN)[avail] == curr.deg
    if(sum(candidates) > 0){
      flag <- T
      not.rec[i] <- ifelse(sum(candidates) > 1, 
                           sample(avail[candidates], 1), 
                           avail[candidates])
    }else{
      curr.deg <- curr.deg + 1
    }
  }
}

degs_tf <- as.data.frame(qqplot(degree(hPIN)[tf], degree(hPIN)[not.tf], 
                                plot.it = F))
degs_rec <- as.data.frame(qqplot(degree(hPIN)[rec], degree(hPIN)[not.rec], 
                                 plot.it = F))
p_tf <- ggplot(degs_tf, aes(x, y)) + geom_point() +
  labs(x = "Degree of TFs", y = "Degree of non-TFs") + theme_bw()
p_rec <- ggplot(degs_rec, aes(x, y)) + geom_point() +
  labs(x = "Degree of receptors", y = "Degree of non-receptors") + theme_bw()

save(not.rec, not.tf, p_tf, p_rec, file = "../results/nonTF_nonRec_pool.RData")
