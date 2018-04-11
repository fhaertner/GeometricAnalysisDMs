
library(dplyr)

chunks <- 32

ref_eff_final <- 0
ref_hop_final <- 0
eff_final <- list()
imp_final <- tibble()

for(i in 1:chunks){
  load(paste0("results/nav_impact_", i, ".RData"))
  ref_eff_final <- ref_eff_final + ref_eff
  ref_hop_final <- ref_hop_final + ref_hop
  eff_final <- c(eff_final, eff)
  if(i == 1){
    imp_final <- imp
  }else{
    imp_final <- rbind(imp_final, imp)
  }
}

ref_eff <- ref_eff_final/chunks
ref_hop <- ref_hop_final/chunks
eff <- eff_final
imp <- imp_final


save(ref_eff, ref_hop, eff, imp, file = "results/nav_impact.RData")