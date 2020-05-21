// Stan program to estimate a  JPLP process
// This is a simple case, assuming multiple shifts come from the same driver
// Different from NHPP with PLP intensity function, in which the likelihood function was evaluated by shifts, this JPLP likelihood function is evaluated by trips. In this way, the likelihood function can be evaluated using the `segment` function in Stan.
functions{
  // LogLikelihood function for shifts with events (N_{event} > 0)
  real jplp_log(vector t_event, vector t_stop, real beta, real theta, real tau, real kappa){
    vector[num_elements(t_event)] loglik_P_1;
    vector[num_elements(t_stop)] loglik_P_2;
    real loglikelihood;
    for (i in 1:num_elements(t_event)){
      loglik_P_1[i] = log(beta) - beta*log(theta) + (beta - 1)*log(t_event[i]);
    }
    loglikelihood = sum(loglik_P_1) - (tau/theta)^beta;
    return loglikelihood;
  }
  // LogLikelihood function for shifts with no event (N_{event} = 0)
  real jplpoevent_log(real tau, real beta, real theta){
    real loglikelihood = - (tau/theta)^beta;
    return(loglikelihood);
  }
}
data {
  int<lower=0> N; //total # of obs
  int<lower=0> S; //total # of shifts
  vector<lower=0>[S] tau;//truncated time
  vector<lower=0>[N] event_time; //failure time
  int group_size[S]; //group sizes
}
parameters{
  real<lower=0> beta;
  real<lower=0> theta;
  real<lower=0> kappa;
}
model{
  int position = 1;
  vector[S] theta_temp;

  // This part should be rewritten to be based on trips, not shifts
  for (s1 in 1:S){
    if(group_size[s1] == 0) {
      tau[s1] ~ jplpoevent_log(beta, theta_temp[s1]);
      }else{
      segment(event_time, position, group_size[s1]) ~ jplp_log(beta, theta_temp[s1], tau[s1], kappa);
      position += group_size[s1];
    }
  }
//PRIORS
  beta ~ gamma(1, 1);
  theta ~ gamma(1, 0.01);
}
