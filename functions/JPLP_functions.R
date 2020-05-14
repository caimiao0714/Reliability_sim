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
