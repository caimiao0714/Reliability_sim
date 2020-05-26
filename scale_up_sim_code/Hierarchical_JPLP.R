pacman::p_load(rstan, dplyr, data.table, broom)
source("functions/JPLP_functions.R")
N_sim = 1000

set.seed(123)
dt1 = sim_hier_JPLP(beta = 1.2,
                    kappa = 0.8,
                    mu0 = 0.2,
                    sigma0 = 0.5,
                    R_K = c(1, 0.3, 0.2),
                    group_size_lambda = 50,
                    D = 10)
fit = stan("stan/jplp_hierarchical.stan",
           chains = 1, iter = 3000, refresh = 0,
           data = dt1$stan_dt, seed = 123)
f_result = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit)




# ----------- Start of simulation -----------------
# D = 10
set.seed(123)
sim10 = list()
for (i in 1:N_sim) {
  print(paste0("D = 10: ", i, " (out of 1000)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, D = 10)
  fit0 = stan("stan/jplp_hierarchical.stan",
              chains = 1, iter = 3000, refresh = 0,
              data = z$stan_dt, seed = 123
  )}, error=function(e){})

  sim10[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim10),
                   "fit/JPLP_fit_sim_hierarchical/sim10.csv")

#########################################################
################         D = 25       ###################
#  ------------------D = 25: Part 1----------------------
set.seed(123)
sim25 = list()
for (i in 1:500) {
  print(paste0("D = 25 (p1), progress: ",
               round(i*100/500, 2),
               "% (", i, " out of 500)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8,
                              mu0 = 0.2, sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2), D = 25)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim25[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim25),
                   "fit/JPLP_fit_sim_hierarchical/sim25_p1.csv")

#  ------------------D = 25: Part 1----------------------
set.seed(124)
sim25 = list()
for (i in 1:500) {
  print(paste0("D = 25 (p2), progress: ",
               round(i*100/500, 2),
               "% (", i, " out of 500)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8,
                              mu0 = 0.2, sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2), D = 25)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim25[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim25),
                   "fit/JPLP_fit_sim_hierarchical/sim25_p2.csv")



#########################################################
################         D = 50       ###################
#  ------------------D = 50: Part 1----------------------
set.seed(123)
sim50 = list()
for (i in 1:500) {
  print(paste0("D = 50 (p1), progress: ",
               round(i*100/500, 2),
               "% (", i, " out of 500)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8,
                              mu0 = 0.2, sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2), D = 50)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim50[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim50),
                   "fit/JPLP_fit_sim_hierarchical/sim50_p1.csv")

#  ------------------D = 50: Part 2----------------------
set.seed(123)
sim50 = list()
for (i in 1:500) {
  print(paste0("D = 50 (p2), progress: ",
               round(i*100/500, 2),
               "% (", i, " out of 500)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8,
                              mu0 = 0.2, sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2), D = 50)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim50[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim50),
                   "fit/JPLP_fit_sim_hierarchical/sim50_p2.csv")



#########################################################
################         D = 75       ###################
#  ------------------D = 75: Part 1----------------------
set.seed(123)
sim75 = list()
for (i in 1:300) {
  print(paste0("D = 75 (p1), progress: ",
               round(i*100/300, 2),
               "% (", i, " out of 300)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8,
                              mu0 = 0.2, sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2), D = 75)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim75[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim75),
                   "fit/JPLP_fit_sim_hierarchical/sim75_p1.csv")
#  ------------------D = 75: Part 2----------------------
set.seed(124)
sim75 = list()
for (i in 1:400) {
  print(paste0("D = 75 (p2), progress: ", round(i*100/400, 2), "% (", i, " out of 300)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8,
                              mu0 = 0.2, sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2), D = 75)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim75[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim75),
                   "fit/JPLP_fit_sim_hierarchical/sim75_p2.csv")
#  ------------------D = 75: Part 3----------------------
set.seed(125)
sim75 = list()
for (i in 1:400) {
  print(paste0("D = 75 (p3), progress: ", round(i*100/400, 2), "% (", i, " out of 300)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8,
                              mu0 = 0.2, sigma0 = 0.5,
                              R_K = c(1, 0.3, 0.2), D = 75)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim75[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim75),
                   "fit/JPLP_fit_sim_hierarchical/sim75_p3.csv")


#########################################################
################         D = 100     ####################
#  ------------------D = 100: Part 1---------------------
set.seed(123)
sim100 = list()
for (i in 1:300) {
  print(paste0("D = 100 (p1), progress: ", round(i*100/300, 2), "% (", i, " out of 300)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8, mu0 = 0.2,
                              sigma0 = 0.5, R_K = c(1, 0.3, 0.2), D = 100)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim100[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim100),
                   "fit/JPLP_fit_sim_hierarchical/sim100_p1.csv")

#  ------------------D = 100: Part 2---------------------
set.seed(124)
sim100 = list()
for (i in 1:300) {
  print(paste0("D = 100 (p2), progress: ", round(i*100/300, 2), "% (", i, " out of 300)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8, mu0 = 0.2,
                              sigma0 = 0.5, R_K = c(1, 0.3, 0.2), D = 100)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim100[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim100),
                   "fit/JPLP_fit_sim_hierarchical/sim100_p2.csv")

#  ------------------D = 100: Part 3---------------------
set.seed(125)
sim100 = list()
for (i in 1:300) {
  print(paste0("D = 100 (p3), progress: ", round(i*100/300, 2), "% (", i, " out of 300)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8, mu0 = 0.2,
                              sigma0 = 0.5, R_K = c(1, 0.3, 0.2), D = 100)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim100[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim100),
                   "fit/JPLP_fit_sim_hierarchical/sim100_p3.csv")

#  ------------------D = 100: Part 4---------------------
set.seed(126)
sim100 = list()
for (i in 1:300) {
  print(paste0("D = 100 (p4), progress: ", round(i*100/300, 2), "% (", i, " out of 300)"))

  tryCatch({z = sim_hier_JPLP(beta = 1.2, kappa = 0.8, mu0 = 0.2,
                              sigma0 = 0.5, R_K = c(1, 0.3, 0.2), D = 100)},
           error=function(e){})

  tryCatch({fit0 = stan("stan/jplp_hierarchical.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123)},
           error=function(e){})

  sim100[[i]] = pull_use("beta|kappa|mu0_true|sigma0|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim100),
                   "fit/JPLP_fit_sim_hierarchical/sim100_p4.csv")




