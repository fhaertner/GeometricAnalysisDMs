
# Given an igraph network, its temperature and coordinates, evaluates the
# embedding to H2

library(NetHypGeom)
library(cowplot)
library(dplyr)

# Load network and coordinates --------------------------------------------
outname <- "emb_eval_hpin"
#net <- readRDS("data/hint.rds")
#load("results/coords_hPIN_150k.RData")
#coords_lh <- coords

load("data/hPIN_150k.RData")
load("data/coords_hPIN_150k.RData")
net <- hPIN
coords_lh <- coords


# Network properties ------------------------------------------------------

N <- vcount(net)
avg_k <- mean(degree(net))
gma <- fit_power_law(degree(net))$alpha
beta <- 1 / (gma - 1) # Parameter controlling popularity fading
m <- round(avg_k/2) # Parameter controlling average node degree

# Connection probability --------------------------------------------------
conn <- get_conn_probs(net, coords_lh, bins = 20)
theo <- get_theoretical_conn_probs(conn$dist, N, avg_k, gma, Temp)

conn <- rbind(conn, theo)
conn$label <- rep(c("LaBNE+HM", "Theory"), each = 20)

p_conn <- ggplot(conn, aes(dist, prob+0.00001, colour = label, shape = label)) + 
  geom_point(size = 3) + geom_line() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(sides = "l") +
  scale_colour_manual(values = c("#339e2b", "#e3196a")) + 
  labs(x = "Hyperbolic distance", y = "Connection probability") + 
  theme_bw() + theme(legend.title = element_blank(), 
                     legend.background = element_blank(), 
                     legend.justification = c(0, 0), legend.position = c(0, 0))


# Real degrees vs expected degrees ----------------------------------------
degs <- tibble(k = degree(net), exp_k = numeric(N))

R <- 2*log(N) - 
  2*log((2*Temp*(1 - exp(-0.5*(1 - beta)*2*log(N))))/(sin(Temp*pi)*m*(1 - beta)))

for(i in 1:N){
  # Compute the hyperbolic distance between a node and all others
  d <- hyperbolic_dist(coords_lh[i, ], coords_lh)
  # Compute the probability that the node will connect with all others
  prob <- 1 / (1 + exp((d - R)/(2*Temp)))
  # Compute the expected node degree
  degs$exp_k[i] <- sum(prob)
}

p_deg <- ggplot(degs, aes(k, round(exp_k))) + geom_point(size = 0.2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format())) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format())) +
  annotation_logticks() + labs(x = "Node degree", y = "Expected node degree") + 
  theme_bw()


# Clustering and routeing -------------------------------------------------

epochs <- 3

cc <- transitivity(net, "average")
theo_cc <- numeric(epochs)

# Source and target nodes based on linear indexing
st <- 1000
idx <- sample(N * N, st)
src <- ((idx - 1) %% N) + 1
trg <- floor((idx - 1) / N) + 1

stretches <- greedy_route_packets(net, coords_lh, src, trg) 
gr <- sum(stretches > 0)/st
hs <- mean(stretches[stretches > 0])
theo_gr <- numeric(epochs)
theo_hs <- numeric(epochs)

for(i in 1:epochs){
  ps_net <- ps_model(N = N, avg.k = avg_k, gma = gma, Temp = Temp)
  theo_cc[i] <- transitivity(ps_net$network, "average")
  stretches <- greedy_route_packets(ps_net$network, ps_net$polar, src, trg)
  theo_gr[i] <- sum(stretches > 0)/st
  theo_hs[i] <- mean(stretches[stretches > 0])
}

clust <- tibble(label = c("Real", "Theory"), cc = c(cc, mean(theo_cc)),
                err = c(0, sd(theo_cc)))

lbl <- factor(c("Greedy routing (GR)", "GR Theory",
                "Hop stretch (HS)", "HS Theory"), 
              levels = c("Greedy routing (GR)", "GR Theory",
                         "Hop stretch (HS)", "HS Theory"),
              ordered = TRUE)
routeing <- tibble(label = lbl, 
                   gr = c(gr, mean(theo_gr), hs, mean(theo_hs)), 
                   err = c(0, sd(theo_gr), 0, sd(theo_hs)))

dodge <- position_dodge(width = 0.9)
p_cc <- ggplot(clust, aes(label, cc)) + geom_col(position = dodge, width = 0.5) +
  geom_errorbar(aes(ymin = cc - err, ymax = cc + err), 
                position = dodge, width = 0.25) + 
  labs(x = "", y = "Clustering coefficient") + theme_bw()

p_gr <- ggplot(routeing, aes(label, gr)) + geom_col(position = dodge, width = 0.5) +
  geom_errorbar(aes(ymin = gr - err, ymax = gr + err), 
                position = dodge, width = 0.25) + 
  labs(x = "", y = "Success rate (%) / Hop stretch") + theme_bw()

fig <- plot_grid(p_conn, p_deg, p_cc, p_gr, nrow = 2, ncol = 2, labels = letters[1:4])

save(fig, conn, degs, clust, routeing,
     file = paste0("results/", outname, ".RData"))

