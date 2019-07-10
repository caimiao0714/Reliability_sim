# simulating PLP - time truncated case
sim_plp_tau = function(tau = 30,
                       beta = 2,
                       theta = 10){
  # initialization
  s = 0; t = 0
  while (max(t) <= tau) {
    u <- runif(1)
    s <- s - log(u)
    t_new <- theta*s^(1/beta)
    t <- c(t, t_new)
  }
  t = t[c(-1, -length(t))]

  return(t)
}

# simulating PLP - failure truncated case
sim_plp_n = function(mean_n, beta, theta){
  N = rpois(1, mean_n)
  u = runif(N, 0, 1)
  n_logu = -log(u)
  s = cumsum(n_logu)
  Delta_t = theta*s^(1/beta)
  return(Delta_t)
}

# simulate multiple NHPPs - failure truncated case
sim_mul_plp_n = function(n_shift = 20,
                         shift_len_mean = 20, shift_len_sd = 5,
                         theta = 10, beta = 2, mean_n = 5){
  #end_time = rnorm(n_shift, shift_len_mean, shift_len_sd)#to be deleted
  t_list = list()
  len_list = list()
  end_time1 = list()

  for (i in 1:n_shift) {
    t_list[[i]] = sim_plp_n(mean_n, beta, theta)
    len_list[[i]] = length(t_list[[i]])
    end_time1[[i]] = ifelse(length(t_list[[i]]) == 0, 0,
                            max(t_list[[i]]))
  }

  event_dat = data.frame(
    shift_id = rep(1:n_shift, unlist(len_list)),
    event_time = Reduce(c, t_list)
  )

  start_end_dat = data.frame(
    shift_id = 1:n_shift,
    start_time = rep(0, n_shift),
    end_time = Reduce(c, end_time1)#to be changed
  )
  start_end_dat$end_time[start_end_dat$end_time == 0] =
    mean(start_end_dat$end_time) # deal with no events

  return(list(event_dat = event_dat,
              start_end_dat = start_end_dat,
              shift_length = unlist(len_list)))
}

# simulate multiple NHPPs - time truncated case
sim_mul_plp_tau = function(n_shift = 20,
                           shift_len_mean = 20, shift_len_sd = 5,
                           theta = 10, beta = 2, mean_n = 5){
  tau_vector = rnorm(n_shift, shift_len_mean, shift_len_sd)#difference1

  t_list = list()
  len_list = list()
  # end_time1 = list() # not needed for time truncated case


  for (i in 1:n_shift) {
    t_list[[i]] = sim_plp_tau(tau_vector[i], beta, theta)
    len_list[[i]] = length(t_list[[i]])
  }

  event_dat = data.frame(
    shift_id = rep(1:n_shift, unlist(len_list)),
    event_time = Reduce(c, t_list)
  )

  start_end_dat = data.frame(
    shift_id = 1:n_shift,
    start_time = rep(0, n_shift),
    end_time = tau_vector #difference2
  )

  return(list(event_dat = event_dat,
              start_end_dat = start_end_dat,
              shift_length = unlist(len_list)))
}

plot_events = function(event_dat, start_end_dat, cross_size = 2){
  p = event_dat %>%
    ggplot(aes(x = event_time, y = shift_id)) +
    geom_point(alpha = 0.8, shape = 4, color = 'red', size = cross_size) +
    scale_y_continuous("shift ID",
                       labels = as.character(start_end_dat$shift_id),
                       breaks = start_end_dat$shift_id)+
    xlab('Time to event (minutes)') +
    geom_segment(data = start_end_dat,
                 aes(x = start_time, xend = end_time,
                     y = shift_id, yend = shift_id),
                 lineend = 'butt',
                 arrow = arrow(length = unit(0.2, "cm"))) +
    theme_classic()
  return(p)
}


sim_hier_plp_tau = function(){
  t_list = list()
  len_list = list()
  tau_vector = rnorm(N, 10, 1.3)

  for (i in 1:N) {
    t_list[[i]] = sim_plp_tau(tau_vector[i], beta, theta[i])
    len_list[[i]] = length(t_list[[i]])
  }

  event_dat = data.frame(
    shift_id = rep(1:N, unlist(len_list)),
    event_time = Reduce(c, t_list)
  )

  start_end_dat = data.frame(
    shift_id = 1:N,
    start_time = rep(0, N),
    end_time = tau_vector #difference2
  )

  return(list(event_dat = event_dat,
              start_end_dat = start_end_dat,
              shift_length = unlist(len_list)))
}