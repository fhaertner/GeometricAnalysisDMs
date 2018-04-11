
library(NetHypGeom)
library(dplyr)

load("../data/hPIN_150k.RData")
load("../data/coords_hPIN_150k.RData")
load("../data/protein_classes.RData")

epochs <- 100 #Number of experiments to run
st <- 500 #Number of source-target pairs

N <- vcount(hPIN) #Number of proteins in the network

#All TFs and Receptors
all.tf <- sort(which(V(hPIN)$symbol %in% protein.classes$symbol[protein.classes$tf]))
all.rec <- sort(which(V(hPIN)$symbol %in% protein.classes$symbol[protein.classes$rec]))
overlap <- intersect(all.rec, all.tf)
all.rec <- c(all.rec[!(all.rec %in% overlap)], overlap)
all.tf <- c(overlap, all.tf[!(all.tf %in% overlap)])

#To form a non-redundant pool of TF-Rec pairs, we use zero-based row-order indexing of a rectangular matrix
n <- length(all.rec)
m <- length(all.tf)

pool.real <- 0:(n*m - 1)
rec.matrix <- matrix(pool.real, byrow = T, nrow = n, ncol = m, dimnames = list(0:(n-1), 0:(m-1)))
overlap.size <- length(overlap)

unwanted <- rec.matrix[as.character((n-overlap.size):(n-1)), as.character(0:(overlap.size-1))]
unwanted <- unwanted[lower.tri(unwanted, diag = T)]
pool.real <- pool.real[-(unwanted+1)]

#Proteins that are neither TFs nor Receptors but with degree similar to them
load("../results/nonTF_nonRec_pool.RData")
not.rec <- sort(not.rec)
not.tf <- sort(not.tf)
overlap <- intersect(not.rec, not.tf)
not.rec <- c(not.rec[!(not.rec %in% overlap)], overlap)
not.tf <- c(overlap, not.tf[!(not.tf %in% overlap)])

nn <- length(not.rec)
mm <- length(not.tf)

pool.ctl <- 0:(nn*mm - 1)
rec.matrix <- matrix(pool.ctl, byrow = T, nrow = nn, ncol = mm, dimnames = list(0:(nn-1), 0:(mm-1)))
overlap.size <- length(overlap)

unwanted <- rec.matrix[as.character((nn-overlap.size):(nn-1)), as.character(0:(overlap.size-1))]
unwanted <- unwanted[lower.tri(unwanted, diag = T)]
pool.ctl <- pool.ctl[-(unwanted+1)]

rm(rec.matrix)


hops.h2 <- vector("list", length = epochs) #Packet delivery using hyperbolic distances
hops.rt <- vector("list", length = epochs) #Packet delivery between Receptors and TF
hops.nrt <- vector("list", length = epochs) #Packet delivery between proteins with degrees similar to Rec and TF

hops.h2 <- vector("list", length = epochs) #Packet delivery using hyperbolic distances
hops.rt <- vector("list", length = epochs) #Packet delivery between Receptors and TF
hops.nrt <- vector("list", length = epochs) #Packet delivery between proteins with degrees similar to Rec and TF

for(ep in 1:epochs){
  #Sampling of non-redundant src-trg pairs
  k <- sample(N*(N-1)/2, st) - 1 #We subtract 1, because the formulae to go from linear upper diagonal indexing to (i,j) are zero-based
  sources <- N - 2 - floor(sqrt(-8*k + 4*N*(N-1)-7)/2.0 - 0.5)
  targets <- k + sources + 1 - N*(N-1)/2 + (N-sources)*((N-sources)-1)/2
  #We sum 1 to go back to 1-based indexing
  sources <- sources + 1
  targets <- targets + 1
  
  #Sampling of non.redundant Rec-TF pairs
  k <- sample(pool.real, st)
  i <- floor(k/m)
  j <- k - (i*m)
  rec <- all.rec[i + 1]
  tf <- all.tf[j + 1]
  
  #Sampling of non.redundant nonRec-nonTF pairs
  k <- sample(pool.ctl, st)
  i <- floor(k/mm)
  j <- k - (i*mm)
  nrec <- not.rec[i + 1]
  ntf <- not.tf[j + 1]
  
  hops.h2[[ep]] <- greedy_route_packets(hPIN, coords, sources, targets)

  hops.rt[[ep]] <- greedy_route_packets(hPIN, coords, rec, tf)
 
  hops.nrt[[ep]] <- greedy_route_packets(hPIN, coords, nrec, ntf)
}

res <- tibble(case = factor(rep(c("Rnd. src-trg pairs", "Rec-TF", "Control"), 
                                each = epochs*2), 
                            levels = c("Rnd. src-trg pairs", "Rec-TF", "Control"), 
                            ordered = TRUE), 
              exp = rep(rep(c("GR efficiency", "Hop stretch"), each = epochs), 3), 
              value = c(sapply(hops.h2, function(x) sum(x > 0)/st),
                        sapply(hops.h2, function(x) mean(x[x > 0])),
                        sapply(hops.rt, function(x) sum(x > 0)/st),
                        sapply(hops.rt, function(x) mean(x[x > 0])),
                        sapply(hops.nrt, function(x) sum(x > 0)/st),
                        sapply(hops.nrt, function(x) mean(x[x > 0]))))

rt <- filter(res, case == "Rec-TF" & exp == "GR efficiency")$value
rnd <- filter(res, case == "Rnd. src-trg pairs" & exp == "GR efficiency")$value
ctrl <- filter(res, case == "Control" & exp == "GR efficiency")$value
pval_RTvsRnd_eff <- wilcox.test(rt, rnd, alternative = "greater")$p.value
pval_RTvsCtl_eff <- wilcox.test(rt, ctrl, alternative = "greater")$p.value

rt <- filter(res, case == "Rec-TF" & exp == "Hop stretch")$value
rnd <- filter(res, case == "Rnd. src-trg pairs" & exp == "Hop stretch")$value
ctrl <- filter(res, case == "Control" & exp == "Hop stretch")$value
pval_RTvsRnd_hs <- wilcox.test(rt, rnd, alternative = "less")$p.value
pval_RTvsCtl_hs <- wilcox.test(rt, ctrl, alternative = "less")$p.value


p.nav <- ggplot(res, aes(case, value, fill = exp)) + geom_boxplot(width = 0.7) +
  labs(x = "", y = "") + theme_bw() +
  theme(legend.title = element_blank(), legend.background = element_blank(), 
        legend.justification = c(1,0.5), legend.position = c(1,0.5))

save(hops.h2, hops.rt, hops.nrt, res, p.nav, pval_RTvsCtl_eff, pval_RTvsCtl_hs, 
     pval_RTvsRnd_eff, pval_RTvsRnd_hs, file = "../results/nav_hPIN.RData")
