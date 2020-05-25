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

pull_use = function(var = "theta", est_obj = f){
  z = est_obj %>%
    broom::tidy() %>%
    filter(grepl(var, term))
  return(z)
}