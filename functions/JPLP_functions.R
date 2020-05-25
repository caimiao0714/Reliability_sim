# sample the number of stops from 1:4
get_n_stop = function() sample(1:4, 1, TRUE)

# Define a inverse function for mean function Lambda
inverse = function (f, lower = 0.0001, upper = 1000) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}


# Mean function Lambda for PLP
Lambda_PLP = function(t, beta = 1.5, theta = 4) return((t/theta)^beta)

# Mean function Lambda for JPLP
Lambda_JPLP = function(t,
                       tau = 12,
                       kappa = 0.8,
                       t_trip = c(3.5, 6.2, 9),
                       beta = 1.5,
                       theta = 4)
{
  t_trip1 = c(0, t_trip)
  n_trip = length(t_trip1)
  comp = Lambda_PLP(t_trip, beta, theta)
  kappa_vec0 = rep(kappa, n_trip - 1)^(0:(n_trip - 2))
  kappa_vec1 = rep(kappa, n_trip - 1)^(1:(n_trip - 1))
  cum_comp0 = comp*kappa_vec0
  cum_comp1 = comp*kappa_vec1
  index_trip = max(cumsum(t > t_trip1)) - 1

  if(index_trip == 0){
    return((t/theta)^beta)
  }else{
    return(sum(cum_comp0[1:index_trip]) - sum(cum_comp1[1:index_trip]) +
             kappa^index_trip*(t/theta)^beta)
  }
}

# Simulate for a JPLP
sim_jplp = function(tau0 = 12,
                    kappa0 = 0.8,
                    t_trip0 = c(3.5, 6.2, 9),
                    beta0 = 1.2,
                    theta0 = 0.5)
{ s = 0; t = 0
  Lambda1 = function(t, tau1 = tau0, kappa1 = kappa0, t_trip1 = t_trip0,
                     beta1 = beta0, theta1 = theta0)
    {
    return(Lambda_JPLP(t, tau = tau1, kappa = kappa1, t_trip = t_trip1,
                  beta = beta1, theta = theta1))
    }
  inv_Lambda = inverse(Lambda1, 0.0001, 1000)

  while (max(t) <= tau0)
    {
    u <- runif(1)
    s <- s - log(u)
    t_new <- inv_Lambda(s)$root
    t <- c(t, t_new)
    }

  t = t[c(-1, -length(t))]

  return(t)
}



# Simulate event times for multiple shifts
sim_mul_jplp = function(kappa = 0.8, beta = 1.2, theta = 2, n_shift = 10)
{
  t_shift_vec = list()
  n_trip_vec = list()
  id_trip_vec = list()
  t_start_vec = list()
  t_stop_vec = list()
  n_event_shift_vec = list()
  n_event_trip_vec = list()
  t_event_vec = list()

  for (i in 1:n_shift) {
    sim_tau = rnorm(1, 10, 1.3)
    n_stop = get_n_stop()
    sim_t_trip = round((1:n_stop)*sim_tau/(n_stop + 1) +
                         rnorm(n_stop, 0, sim_tau*0.15/n_stop), 2)
    t_events = sim_jplp(tau0 = sim_tau,
                        kappa0 = kappa,
                        t_trip0 = sim_t_trip,
                        beta0 = beta,
                        theta0 = theta)
    t_shift_vec[[i]] = sim_tau
    id_trip_vec[[i]] = 1:(n_stop + 1)
    n_trip_vec[[i]] = n_stop + 1
    t_start_vec[[i]] = c(0, sim_t_trip)
    t_stop_vec[[i]]  = c(sim_t_trip, sim_tau)
    n_event_shift_vec[[i]] = length(t_events)
    t_event_vec[[i]] = t_events

    tmp_n_event_trip = rep(NA_integer_, (n_stop + 1))
    for (j in 1:(n_stop + 1)) {
       tmp_n_event_trip[j] = sum(t_events > t_start_vec[[i]][j] &
                                   t_events <= t_stop_vec[[i]][j])
    }
    n_event_trip_vec[[i]] = tmp_n_event_trip
  }

  event_dt = data.frame(
    shift_id = rep(1:n_shift, unlist(n_event_shift_vec)),
    trip_id = rep(Reduce(c, id_trip_vec), Reduce(c, n_event_trip_vec)),
    event_time = Reduce(c, t_event_vec)
  )

  trip_dt = data.frame(
    shift_id = rep(1:n_shift, Reduce(c, n_trip_vec)),
    trip_id = Reduce(c, id_trip_vec),
    t_trip_start = Reduce(c, t_start_vec),
    t_trip_end = Reduce(c, t_stop_vec),
    N_events = Reduce(c, n_event_trip_vec)
  )

  shift_dt = data.frame(
    shift_id = 1:n_shift,
    start_time = rep(0, n_shift),
    end_time = Reduce(c, t_shift_vec)
  )

  stan_dt = list(N = nrow(event_dt),
                 S = nrow(trip_dt),
                 r_trip = trip_dt$trip_id,
                 t_trip_start = trip_dt$t_trip_start,
                 t_trip_end = trip_dt$t_trip_end,
                 event_time = event_dt$event_time,
                 group_size = trip_dt$N_events)

  return(list(event_dt = event_dt,
              trip_dt = trip_dt,
              shift_dt = shift_dt,
              stan_dt = stan_dt))

}

sim_hier_JPLP = function(
  beta = 1.2,
  kappa = 0.8,
  D = 10, # the number of drivers
  K = 3, # the number of predictor variables
  group_size_lambda = 10, # the mean number of shifts for each driver
  mu0 = 0.2, # hyperparameter 1
  sigma0 = 0.5, # hyperparameter 2
  R_K = c(1, 0.3, 0.2) # Fixed-effects parameters
)
{
  # 1. Random-effect intercepts
  r_0D = rnorm(D, mean = mu0, sd = sigma0)

  # 3. The number of observations (shifts) in the $d$-th driver: $N_{d}$
  N_K = rpois(D, group_size_lambda)
  N = sum(N_K) # the total number of shifts for all D drivers
  id = rep(1:D, N_K)

  # 4. Generate data: x_1, x_2, .. x_K
  simX = function(group_sizes = N_K)
  {
    ntot = sum(group_sizes)

    int1 = rep(1, ntot)
    x1 =  rnorm(ntot, 1, 1)
    x2 = rgamma(ntot, 1, 1)
    x3 =  rpois(ntot, 2)

    return(data.frame(int1, x1, x2, x3))
  }
  X = simX(N_K)

  # 5. Scale parameters of a JPLP
  # 5a. parameter matrix: P
  P = cbind(r0 = rep(r_0D, N_K), t(replicate(N, R_K)))
  M_logtheta = P*X
  theta = exp(rowSums(M_logtheta))

  # Initialization of lists
  t_shift_vec = list()
  n_stop_vec = list()
  t_stop_vec = list()
  n_event_vec = list()
  t_event_vec = list()

  for (i in 1:N)
  {
    sim_tau = rnorm(1, 10, 1.3)
    n_stop = get_n_stop()
    sim_t_trip = round((1:n_stop)*sim_tau/(n_stop + 1) +
                         rnorm(n_stop, 0, sim_tau*0.15/n_stop), 2)
    t_events = sim_jplp(tau0 = sim_tau,
                        kappa0 = kappa,
                        t_trip0 = sim_t_trip,
                        beta0 = beta,
                        theta0 = theta[i])
    t_shift_vec[[i]] = sim_tau
    n_stop_vec[[i]] = n_stop
    t_stop_vec[[i]] = sim_t_trip
    n_event_vec[[i]] = length(t_events)
    t_event_vec[[i]] = t_events
  }



  # shifts data
  shift_dt = data.frame(
    driver_id = rep(1:D, N_K),
    shift_id = 1:N,
    start_time = rep(0, N),
    end_time = Reduce(c, t_shift_vec),
    n_stop = Reduce(c, n_stop_vec),
    n_event = Reduce(c, n_event_vec)
  )

  # trips data set
  trip_dt = data.frame(
    driver_id = rep(shift_dt$driver_id, shift_dt$n_stop),
    shift_id = rep(1:N, unlist(n_stop_vec)),
    trip_time = Reduce(c, t_stop_vec)
  )

  # TEMPORARY vector: a temporary vector for events per driver
  n_event_driver = shift_dt %>%
    group_by(driver_id) %>%
    summarise(n_event = sum(n_event)) %>%
    pull(n_event)

  # events data set
  event_dt = data.frame(
    driver_id = rep(1:D, n_event_driver),
    shift_id = rep(1:N, Reduce(c, n_event_vec)),
    event_time = Reduce(c, t_event_vec)
  )

  return(list(event_time = event_dt,
              trip_time = trip_dt,
              shift_time = shift_dt))
}

pull_use = function(var = "theta", est_obj = f){
  z = est_obj %>%
    broom::tidy() %>%
    filter(grepl(var, term))
  return(z)
}