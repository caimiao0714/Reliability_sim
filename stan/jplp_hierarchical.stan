// Stan program to estimate a  JPLP process
// This is a simple case, assuming multiple shifts come from the same driver
// Different from NHPP with PLP intensity function, in which the likelihood function was evaluated by shifts, this JPLP likelihood function is evaluated by trips. In this way, the likelihood function can be evaluated using the `segment` function in Stan.
functions{
  // LogLikelihood function for shifts with events (N_{event} > 0)
  real jplp_log(vector t_event, // time of SCEs
                real trip_start,
                real trip_end,
                int r,// trip index
                real beta,
                real theta,
                real kappa)
  {
    vector[num_elements(t_event)] loglik;
    real loglikelihood;
    for (i in 1:num_elements(t_event))
    {
      loglik[i] = (r - 1)*log(kappa) + log(beta) - beta*log(theta) + (beta - 1)*log(t_event[i]);
    }
    loglikelihood = sum(loglik) - kappa^(r - 1)*theta^(-beta)*(trip_end^beta - trip_start^beta);
    return loglikelihood;
  }
  // LogLikelihood function for shifts with no event (N_{event} = 0)
  real jplpoevent_lp(real trip_start,
                     real trip_end,
                     int r,
                     real beta,
                     real theta,
                     real kappa)
  {
    real loglikelihood = - kappa^(r - 1)*theta^(-beta)*(trip_end^beta - trip_start^beta);
    return(loglikelihood);
  }
}
data {
  int<lower=0> N; //total # of events
  int<lower=1> D; //total # of drivers
  int<lower=1> K; //number of predictors
  int<lower=0> S; //total # of trips
  int<lower=1> id[S];//driver index, must be an array
  int r_trip[S];//index of trip $r$
  vector<lower=0>[S] t_trip_start;//trip start time
  vector<lower=0>[S] t_trip_end;//trip end time
  vector<lower=0>[N] event_time; //failure time
  int group_size[S]; //group sizes
  matrix[S, K] X_predictors;//predictor variable matrix
}
transformed data{
  matrix[S, K] X_centered;
  vector[K] X_means;
  for(k0 in 1:K){
    X_means[k0] = mean(X_predictors[, k0]);
    X_centered[,k0] = X_predictors[, k0] - X_means[k0];
  }
}
parameters{
  real<lower=0> beta;
  real<lower=0> theta;
  real<lower=0> kappa;
  vector[K] R1_K; // fixed parameters for K predictors
  vector[D] R0; // random intercept for D drivers
}
model{
  int position = 1;
  vector[S] theta_temp;

  for (s0 in 1:S){
    theta_temp[s0] = exp(R0[id[s0]] + X_centered[s0,]*R1_K);
  }

  for (s1 in 1:S){ // Likelihood estimation for JPLP based on trips, not shifts
    if(group_size[s1] == 0){
      target += jplpoevent_lp(t_trip_start[s1], t_trip_end[s1], r_trip[s1], beta, theta_temp[s1], kappa);
      }else{
      segment(event_time, position, group_size[s1]) ~ jplp_log(t_trip_start[s1], t_trip_end[s1], r_trip[s1], beta, theta_temp[s1], kappa);
      position += group_size[s1];
    }
  }
//PRIORS
  beta ~ gamma(1, 1);
  theta ~ gamma(1, 0.01);
  kappa ~
}
