
library(cowplot)


# FIGURE S1 ---------------------------------------------------------------

load("results/emb_eval_hpin.RData")

save_plot(filename = "figs/fig_s01.svg", plot = fig, nrow = 2, ncol = 2)

# FIGURE S2 ---------------------------------------------------------------

load("results/bio_eval_hpin.RData")

figS2 <- ggdraw() + 
  draw_plot(hpin.age, x = 0, y = 0, width = 1/3, height = 1) + 
  draw_plot(hpin.tht, x = 1/3, y = 0, width = 2/3, height = 1) +
  draw_plot_label(letters[1:2], x = c(0, 1/3), 
                  y = c(1, 1), size = 20)

save_plot(filename = "figs/fig_s02.svg", plot = figS2, nrow = 1, ncol = 3)

# FIGURE S3 ---------------------------------------------------------------

load("results/dm_angles.RData")

save_plot(filename = "figs/fig_s03.svg", plot = p_supp, nrow = 3, ncol = 2)


# FIGURE S4 ---------------------------------------------------------------

load("results/ndr_separation.RData")

partA <- pS + theme(legend.position = "none")
partB <- pJ + theme(legend.position = "none")

figS4 <- ggdraw() + 
  draw_plot(partA, x = 0, y = 0, width = 0.5, height = 1) + 
  draw_plot(partB, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(letters[1:2], x = c(0, 0.5), 
                  y = c(1, 1), size = 15)

save_plot(filename = "figs/fig_s04.svg", plot = figS4, nrow = 1.7, ncol = 2)


# FIGURE S5 ---------------------------------------------------------------

# Dendrogram


# FIGURE S6 ---------------------------------------------------------------

load("results/nonTF_nonRec_pool.RData")

figS6 <- ggdraw() + 
  draw_plot(p_rec, x = 0, y = 0, width = 0.5, height = 1) + 
  draw_plot(p_tf, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(letters[1:2], x = c(0, 0.5), 
                  y = c(1, 1), size = 15)

save_plot(filename = "figs/fig_s06.svg", plot = figS6, nrow = 1, ncol = 2)


# FIGURE S7 ---------------------------------------------------------------

load("results/go_impact.RData")

figS7 <-  ggdraw() + 
  draw_plot(p_imp, x = 0, y = 0.5, width = 1, height = 0.5) + 
  draw_plot(p_nim, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(letters[1:2], x = c(0, 0), 
                  y = c(1, 0.5), size = 15)

save_plot(filename = "figs/fig_s07.svg", plot = figS7, nrow = 3, ncol = 1.7)


# FIGURE S8 ---------------------------------------------------------------

load("results/nav_imp_drugs.RData")

figS8 <-  ggdraw() + 
  draw_plot(p_freq, x = 0, y = 0.5, width = 1, height = 0.5) + 
  draw_plot(p_infreq, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(letters[1:2], x = c(0, 0), 
                  y = c(1, 0.5), size = 15)

save_plot(filename = "figs/fig_s08.svg", plot = figS8, nrow = 3, ncol = 1.7)

# FIGURE S9 ---------------------------------------------------------------

load("results/emb_eval_hint.RData")

save_plot(filename = "figs/fig_sx1.svg", plot = fig, nrow = 2, ncol = 2)

# FIGURE S10 ---------------------------------------------------------------

load("results/bio_eval_hint.RData")

figS10 <- ggdraw() + 
  draw_plot(hint.age, x = 0, y = 0, width = 1/3, height = 1) + 
  draw_plot(hint.tht, x = 1/3, y = 0, width = 2/3, height = 1) +
  draw_plot_label(letters[1:2], x = c(0, 1/3), 
                  y = c(1, 1), size = 20)

save_plot(filename = "figs/fig_s10.svg", plot = figS10, nrow = 1, ncol = 3)


# FIGURE S11 ---------------------------------------------------------------

load("results/S_dm_stats_hint.RData")

partA <- ggplot(dm_stats, aes(dm_size, comp)) + geom_point() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
                labels = scales::trans_format("log10", scales::math_format())) +
  annotation_logticks(sides = "b") +
  labs(x = "Disease module size", y = "Connected components") +
  theme_bw()

partB <- ggplot(dm_stats, aes(lcc/dm_size, -log10(pval_lcc))) + geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "red") +
  geom_vline(xintercept = 1, linetype = 2, colour = "blue") +
  labs(x = "Fraction of disease proteins in LCC", 
       y = expression(-log[10](p-value))) +
  theme_bw()

partC <- ggplot(dm_stats, aes(lcc/dm_size, -log10(pval_sp))) + geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "red") +
  geom_vline(xintercept = 1, linetype = 2, colour = "blue") +
  labs(x = "Fraction of disease proteins in LCC", 
       y = expression(-log[10](p-value))) +
  theme_bw()

partD <- ggplot(dm_stats, aes(lcc/dm_size, -log10(pval_h2))) + geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "red") +
  geom_vline(xintercept = 1, linetype = 2, colour = "blue") +
  labs(x = "Fraction of proteins in LCC", y = expression(-log[10](p-value))) +
  theme_bw()


figS11 <-  ggdraw() + 
  draw_plot(partA, x = 0, y = 0.5, width = 0.5, height = 0.5) + 
  draw_plot(partB, x = 0.5, y = 0.5, width = 0.5, height = 0.5) + 
  draw_plot(partC, x = 0, y = 0, width = 0.5, height = 0.5) + 
  draw_plot(partD, x = 0.5, y = 0, width = 0.5, height = 0.5) + 
  draw_plot_label(letters[1:4], x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.5, 0.5), size = 15)

save_plot(filename = "figs/fig_s11.svg", plot = figS11, ncol = 2, nrow = 2)


# FIGURE S12 ---------------------------------------------------------------

load("results/dm_angles_hint.RData")

save_plot(filename = "figs/fig_s12.svg", plot = p_supp, nrow = 2, ncol = 2)


# FIGURE S13 ----------------------------------------------------------------

load("results/S_ndr_separation_hint.RData")

save_plot(filename = "figs/fig_S13.svg", plot = pH, nrow = 1.7, ncol = 2)


# FIGURE S14 --------------------------------------------------------------

load("results/S_nav_hint_50epochs.RData")
load("results/nav_imp_drugs.RData")

partA <- ggplot(imp, aes(x = imp_eff_dis, y = imp_eff_rnd, 
                         size = -log10(pval_eff))) + geom_point() + 
  scale_size_continuous(breaks = c(seq(1, 12, 3))) +
  geom_vline(xintercept = min(imp$imp_eff_rnd, na.rm = TRUE), linetype = 2, 
             colour = "blue") +
  labs(x = "Impact of DM proteins on navigability", 
       y = "Impact of rnd. proteins on navigability",
       size = expression(-log[10](p-val))) + 
  theme_bw() + 
  theme(legend.background = element_blank(), 
        legend.justification = c(0,1), legend.position = c(0,1))

partB <- p_dt

figS14 <- plot_grid(partA, partB, labels = letters[1:2])
save_plot(figS14, filename = "figs/fig_s14.svg", nrow = 1, ncol = 2)
