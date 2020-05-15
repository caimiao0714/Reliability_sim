// Stan program to estimate a  JPLP process
// This is a simple case, assuming multiple shifts come from the same driver
functions{
  // LogLikelihood function for shifts with events (N_{event} > 0)
  real jplp_log(vector t, real beta, real theta, real tau, real kappa){
    vector[num_elements(t)] loglik_part;
    real loglikelihood;
    for (i in 1:num_elements(t)){
      loglik_part[i] = log(beta) - beta*log(theta) + (beta - 1)*log(t[i]);
    }
    loglikelihood = sum(loglik_part) - (tau/theta)^beta;
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
