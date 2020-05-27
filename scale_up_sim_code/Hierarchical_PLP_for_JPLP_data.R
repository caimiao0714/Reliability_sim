pacman::p_load(rstan, dplyr, data.table)
source("functions/NHPP_functions.R")
source("functions/JPLP_functions.R")
N_sim = 1000

# test code
set.seed(123)
z = sim_hier_JPLP(beta = 1.2,
                  kappa = 0.8,
                  mu0 = 0.2,
                  sigma0 = 0.5,
                  R_K = c(1, 0.3, 0.2),
                  group_size_lambda = 50,
                  D = 10)
fit0 = stan("stan/nhppnoevent_lp.stan",
            chains = 1, iter = 1000, data = z$stan_jplp_dt_for_plp, refresh = 0)
pull_use(var = "beta|kappa|mu0_true|sigma0|R1_K", fit0)



################################################################
################     Start of simulation     ###################
################################################################
################         D = 10       ##########################
set.seed(123)
sim10 = list()
for (i in 1:N_sim) {
  print(paste0("D = 10: ", round(i*100/N_sim, 2), "% (", i, " out of 1000)"))
  tryCatch({z= sim_hier_JPLP(beta = 1.2,
                             kappa = 0.8,
                             mu0 = 0.2,
                             sigma0 = 0.5,
                             R_K = c(1, 0.3, 0.2),
                             group_size_lambda = 10,
                             D = 10)}, error=function(e){})
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$stan_jplp_dt_for_plp, seed = 123
  )}, error=function(e){})

  sim10[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim10),
                   "fit/PLP_sim_for_JPLP_data/sim10.csv")

################################################################
################         D = 25       ##########################
set.seed(123)
sim25 = list()
for (i in 1:N_sim) {
  print(paste0("D = 25: ", round(i*100/N_sim, 2), "% (", i, " out of 1000)"))
  tryCatch({z = sim_hier_JPLP(beta = 1.2,
                              kappa = 0.8,
                              mu0 = 0.2,
                              sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2),
                              group_size_lambda = 10,
                              D = 25)}, error=function(e){})
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$stan_jplp_dt_for_plp, seed = 123
  )}, error=function(e){})

  sim25[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim25),
                   "fit/PLP_sim_for_JPLP_data/sim25.csv")

################################################################
################         D = 50       ##########################
set.seed(123)
sim50 = list()
for (i in 1:N_sim) {
  print(paste0("D = 50: ", i, " (out of 1000)"))
  tryCatch({z = sim_hier_JPLP(beta = 1.2,
                              kappa = 0.8,
                              mu0 = 0.2,
                              sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2),
                              group_size_lambda = 10,
                              D = 50)}, error=function(e){})
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$stan_jplp_dt_for_plp, seed = 123
  )}, error=function(e){})

  sim50[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim50),
                   "fit/PLP_sim_for_JPLP_data/sim50.csv")


################################################################
################         D = 75       ##########################
set.seed(123)
sim75 = list()
for (i in 1:N_sim) {
  print(paste0("D = 75: ", i, " (out of 1000)"))
  tryCatch({z = sim_hier_JPLP(beta = 1.2,
                              kappa = 0.8,
                              mu0 = 0.2,
                              sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2),
                              group_size_lambda = 10,
                              D = 75)}, error=function(e){})
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$stan_jplp_dt_for_plp, seed = 123
  )}, error=function(e){})

  sim75[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim75),
                   "fit/PLP_sim_for_JPLP_data/sim75.csv")

################################################################
################         D = 100       #########################
set.seed(123)
sim100 = list()
for (i in 1:N_sim) {
  print(paste0("D = 100: ", i, " (out of 1000)"))
  tryCatch({z = sim_hier_JPLP(beta = 1.2,
                              kappa = 0.8,
                              mu0 = 0.2,
                              sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2),
                              group_size_lambda = 10,
                              D = 100)}, error=function(e){})
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$stan_jplp_dt_for_plp, seed = 123
  )}, error=function(e){})

  sim100[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim100),
                   "fit/PLP_sim_for_JPLP_data/sim100.csv")