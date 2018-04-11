
library(cowplot)
library(ggridges)
library(dplyr)
library(stringr)


# FIGURE 1 ----------------------------------------------------------------

load("../results/dm_stats.RData")

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

fig1 <- ggdraw() + 
  draw_plot(partA, x = 0.75, y = 0.75, width = 0.25, height = 0.25) + 
  draw_plot(partB, x = 0.75, y = 0.5, width = 0.25, height = 0.25) +
  draw_plot(partC, x = 0.75, y = 0.25, width = 0.25, height = 0.25) +
  draw_plot(partD, x = 0.75, y = 0.00, width = 0.25, height = 0.25) +
  draw_plot_label(letters[1:5], x = c(0, 0.75, 0.75, 0.75, 0.75), 
                  y = c(1, 1, 0.75, 0.50, 0.25), size = 20)

save_plot(filename = "../figs/fig_01.svg", plot = fig1, nrow = 3, ncol = 3)


# FIGURE 2 ----------------------------------------------------------------

load("results/h2_vs_bp.RData")
load("results/h2_vs_cc.RData")
load("results/dm_angles.RData")

fig2 <- ggdraw() + 
  draw_plot(bp.p, x = 0, y = 0.5, width = 1/3, height = 0.5) + 
  draw_plot(cc.p, x = 0, y = 0, width = 1/3, height = 0.5) +
  draw_plot(p_main, x = 1/3, y = 0, width = 2/3, height = 1) +
  draw_plot_label(letters[1:3], x = c(0, 0, 1/3), 
                  y = c(1, 0.5, 1), size = 15)

save_plot(filename = "../figs/fig_02.svg", plot = fig2, nrow = 2, ncol = 2)


# FIGURE 3 ----------------------------------------------------------------

load("results/ndr_separation.RData")

# Info about some interesting cases
dg %>% 
  filter(diseaseType == "Diseases of the nervous system" & mH2 > 3) %>% 
  select(diseaseName, diseaseType)

dg %>% 
  filter(diseaseType != "Neoplasms" & mH2 > 1 & mH2 < 2) %>% 
  select(diseaseName, mH1) %>% 
  arrange(mH1)

dg %>% 
  filter(diseaseType != "Endocrine, nutritional and metabolic diseases" & 
           mH2 > -2 & mH2 < -1) %>% 
  select(diseaseName, mH1) %>% 
  arrange(mH1)

save_plot(filename = "../figs/fig_03.svg", plot = pH, nrow = 1.7, ncol = 2)


# FIGURE 4 ----------------------------------------------------------------

load("results/nav_hPIN.RData")
load("results/nav_impact_50epochs.RData")
load("results/nav_imp_drugs.RData")

partA <- p.nav + scale_fill_manual(values = c("#ef8a62", "#67a9cf"))

partB <- ggplot(imp, aes(x = imp_eff_dis, y = imp_eff_rnd, 
                         size = -log10(pval_eff))) + geom_point() + 
  scale_size_continuous(breaks = c(seq(1, 15, 3))) +
  geom_vline(xintercept = min(imp$imp_eff_rnd), linetype = 2, colour = "blue") +
  coord_cartesian(xlim = c(-0.05, 0)) +
  labs(x = "Impact of DM proteins on navigability", 
       y = "Impact of rnd. proteins on navigability",
       size = expression(-log[10](p-val))) + 
  theme_bw() + 
  theme(legend.background = element_blank(), 
        legend.justification = c(0,1), legend.position = c(0,1))

partC <- p_dt

fig4 <- ggdraw() + 
  draw_plot(partA, x = 0, y = 0.5, width = 0.5, height = 0.5) + 
  draw_plot(partB, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(partC, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(letters[1:3], x = c(0, 0.5, 0), 
                  y = c(1, 1, 0.5), size = 15)

save_plot(filename = "../figs/fig_04.svg", plot = fig4, nrow = 2, ncol = 2)
