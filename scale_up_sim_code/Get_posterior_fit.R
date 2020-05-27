pacman::p_load(data.table, dplyr, tidyr)
comb_fit = function(file_path){
  csv_files = list.files(file_path)
  csv_files = csv_files[file.size(list.files(file_path, full.names = T)) != 0]

  dt_list = list()
  for (i in seq_along(csv_files)) {
    dt_list[[i]] = fread(paste0(file_path, csv_files[i]))
    dt_list[[i]]$D = as.integer(gsub('sim|_p[[:digit:]]\\.csv|\\.csv', '', csv_files[i]))
  }
  dt = rbindlist(dt_list)
  dt[,sim_scenaio := gsub('fit/|/', '', file_path)]
  return(dt)
}

data.table::fwrite(comb_fit('fit/PLP_sim/'), 'fit/PLP_sim.csv')
data.table::fwrite(comb_fit('fit/JPLP_sim/'), 'fit/JPLP_sim.csv')
data.table::fwrite(comb_fit('fit/PLP_sim_for_JPLP/'), 'fit/PLP_sim_for_JPLP.csv')

d1 = fread('fit/PLP_sim.csv')
d2 = fread('fit/JPLP_sim.csv')
d3 = fread('fit/PLP_sim_for_JPLP.csv')
dt = rbindlist(list(d1, d2, d3))


dt %>%
  .[,.(N_sim = .N,
       estimate_mean = mean(estimate),
       estimate_sd = sd(estimate),
       se_mean = mean(std.error)),
    .(term, D, sim_scenaio)] %>%
  .[,delta_true := case_when(
    term == 'beta' ~ abs(estimate_mean - 1.2),
    term == 'kappa' ~ abs(estimate_mean - 0.8),
    term == 'theta' ~ abs(estimate_mean - 2),
    term == 'mu0_true' ~ abs(estimate_mean - 0.2),
    term == 'sigma0' ~ abs(estimate_mean - 0.5),
    term == 'R1_K[1]' ~ abs(estimate_mean - 1),
    term == 'R1_K[2]' ~ abs(estimate_mean - 0.3),
    term == 'R1_K[3]' ~ abs(estimate_mean - 0.2)
  )] %>%
  .[order(term, D),
    .(sim_scenaio, parameter = term, D, N_sim, estimate_mean, delta_true, #estimate_sd,
      se_mean)] %>%
  .[,sim_scenaio := factor(sim_scenaio, levels = c('PLP_sim', 'JPLP_sim', 'PLP_sim_for_JPLP'))] %>%
  pivot_longer(estimate_mean:se_mean, names_to = "estimate", values_to = "value") %>%
  pivot_wider(id_cols = c(sim_scenaio, D, estimate), names_from = parameter, values_from = value) %>%
  arrange(sim_scenaio, estimate, D) %>%
  knitr::kable('latex', digits = 4, booktabs = T, linesep = "",
               align = c('l', 'r', rep('c', 8)))


