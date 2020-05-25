pacman::p_load(rstan, dplyr, data.table, broom)
source("functions/JPLP_functions.R")
N_sim = 1000


dt = sim_mul_jplp(kappa = 0.8, beta = 1.2, theta = 2, n_shift = 10)
fit = stan("stan/jplp_simple.stan",
           chains = 1, iter = 3000, refresh = 0,
           data = dt$stan_dt, seed = 123)
f_result = pull_use("beta|theta|kappa", fit)


# ----------- Start of simulation -----------------
# n_shift = 10
set.seed(123)
sim10 = list()
for (i in 1:N_sim) {
  print(paste0("N shift = 10: ", i, " (out of 1000)"))
  z = sim_mul_jplp(beta = 1.2, n_shift = 10)
  tryCatch({fit0 = stan("stan/jplp_simple.stan",
                        chains = 1, iter = 3000, refresh = 0,
                        data = z$stan_dt, seed = 123
  )}, error=function(e){})

  sim10[[i]] = pull_use("beta|theta|kappa", fit0)
}
data.table::fwrite(data.table::rbindlist(sim10),
                   "fit/JPLP_sim_simple/sim10.csv")

# n_shift = 25
set.seed(123)
sim25 = list()
for (i in 1:N_sim) {
  print(paste0("N shift = 25: ", i, " (out of 1000)"))
  z = sim_mul_jplp(beta = 1.2, n_shift = 25)
  tryCatch({fit0 = stan("stan/jplp_simple.stan",
                        chains = 1, iter = 3000, refresh = 0,
                        data = z$stan_dt, seed = 123
  )}, error=function(e){})

  sim25[[i]] = pull_use("beta|theta|kappa", fit0)
}
data.table::fwrite(data.table::rbindlist(sim25),
                   "fit/JPLP_sim_simple/sim25.csv")

# n_shift = 50
set.seed(123)
sim50 = list()
for (i in 1:N_sim) {
  print(paste0("N shift = 50: ", i, " (out of 1000)"))

  tryCatch({z = sim_mul_jplp(beta = 1.2, n_shift = 50)
  fit0 = stan("stan/jplp_simple.stan",
              chains = 1, iter = 4000, refresh = 0,
              data = z$stan_dt, seed = 123
  )}, error=function(e){})

  sim50[[i]] = pull_use("beta|theta|kappa", fit0)
}
data.table::fwrite(data.table::rbindlist(sim50),
                   "fit/JPLP_sim_simple/sim50.csv")


# n_shift = 75
set.seed(123)
sim75 = list()
for (i in 1:N_sim) {
  print(paste0("N shift = 75: ", i, " (out of 1000)"))
  z = sim_mul_jplp(beta = 1.2, n_shift = 75)
  tryCatch({fit0 = stan("stan/jplp_simple.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123
  )}, error=function(e){})

  sim75[[i]] = pull_use("beta|theta|kappa", fit0)
}
data.table::fwrite(data.table::rbindlist(sim75),
                   "fit/JPLP_sim_simple/sim75.csv")


# n_shift = 100
set.seed(123)
sim100 = list()
for (i in 1:N_sim) {
  print(paste0("N shift = 100: ", i, " (out of 1000)"))
  z = sim_mul_jplp(beta = 1.2, n_shift = 100)
  tryCatch({fit0 = stan("stan/jplp_simple.stan",
                        chains = 1, iter = 4000, refresh = 0,
                        data = z$stan_dt, seed = 123
  )}, error=function(e){})

  sim100[[i]] = pull_use("beta|theta|kappa", fit0)
}
data.table::fwrite(data.table::rbindlist(sim100),
                   "fit/JPLP_sim_simple/sim100.csv")



# ------------------ Validate estimates ------------
d = fread('fit/JPLP_sim_simple/sim100.csv')
csv_files = list.files('fit/JPLP_sim_simple/')
dt_list = list()
for (i in seq_along(csv_files)) {
  dt_list[[i]] = fread(paste0('fit/JPLP_sim_simple/', csv_files[i]))
  dt_list[[i]]$n_shift = as.integer(gsub('sim|\\.csv', '', csv_files[i]))
}
dt = rbindlist(dt_list)

dt %>%
  .[,.(N = .N,
       estimate = mean(estimate),
       estimate_sd = sd(estimate),
       mean_sd = mean(std.error)),
    .(term, n_shift)] %>%
  .[order(term, n_shift)]












