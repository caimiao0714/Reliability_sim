pacman::p_load(rstan, dplyr, data.table)
source("functions/NHPP_functions.R")
N_sim = 1000


# test code
df = sim_hier_nhpp(D = 10, beta = 1.2)
fit0 = stan("stan/nhppnoevent_lp.stan",
            chains = 1, iter = 1000, data = df$hier_dat, refresh = 0)
f_result = broom::tidy(pull_use)


# --- Start of simulation ---

set.seed(123)
sim10 = list()
for (i in 1:N_sim) {
  print(paste0("D = 10: ", i, " (out of 1000)"))
  z = sim_hier_nhpp(D = 10, beta = 1.2)
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$hier_dat, seed = 123
  )}, error=function(e){})

  sim10[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim10),
                   "fit/NHPP_sim/sim10.csv")


set.seed(123)
sim25 = list()
for (i in 1:N_sim) {
  print(paste0("D = 25: ", i, " (out of 1000)"))
  z = sim_hier_nhpp(D = 25, beta = 1.2)
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$hier_dat, seed = 123
  )}, error=function(e){})

  sim25[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim25),
                   "fit/NHPP_sim/sim25.csv")

set.seed(123)
sim50 = list()
for (i in 1:N_sim) {
  print(paste0("D = 50: ", i, " (out of 1000)"))
  z = sim_hier_nhpp(D = 50, beta = 1.2)
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$hier_dat, seed = 123
  )}, error=function(e){})

  sim50[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim50),
                   "fit/NHPP_sim/sim50.csv")


set.seed(123)
sim75 = list()
for (i in 1:N_sim) {
  print(paste0("D = 75: ", i, " (out of 1000)"))
  z = sim_hier_nhpp(D = 75, beta = 1.2)
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$hier_dat, seed = 123
  )}, error=function(e){})

  sim75[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim75),
                   "fit/NHPP_sim/sim75.csv")


set.seed(123)
sim100 = list()
for (i in 1:N_sim) {
  print(paste0("D = 100: ", i, " (out of 1000)"))
  z = sim_hier_nhpp(D = 100, beta = 1.2)
  tryCatch({fit0 = stan("stan/nhppnoevent_lp.stan",
                        chains = 1, iter = 1000, refresh = 0,
                        data = z$hier_dat, seed = 123
  )}, error=function(e){})

  sim100[[i]] = pull_use("mu0_true|sigma0|beta|R1_K", fit0)
}
data.table::fwrite(data.table::rbindlist(sim100),
                   "fit/NHPP_sim/sim100.csv")


