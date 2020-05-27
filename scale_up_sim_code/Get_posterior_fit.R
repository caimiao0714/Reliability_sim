file_path = 'fit/JPLP_fit_sim_hierarchical/'
csv_files = list.files(file_path)
csv_files = csv_files[file.size(list.files(file_path, full.names = T)) != 0]

dt_list = list()
for (i in seq_along(csv_files)) {
  dt_list[[i]] = fread(paste0('fit/JPLP_fit_sim_hierarchical/', csv_files[i]))
  dt_list[[i]]$D = as.integer(gsub('sim|_p[[:digit:]]\\.csv|\\.csv', '', csv_files[i]))
}
dt = rbindlist(dt_list)

dt %>%
  .[,.(N_sim = .N,
       estimate_mean = mean(estimate),
       estimate_sd = sd(estimate),
       se_mean = mean(std.error)),
    .(term, D)] %>%
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
    .(parameter = term, D, N_sim, estimate_mean, delta_true, 
      estimate_sd, se_mean)]