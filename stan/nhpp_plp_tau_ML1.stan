functions{
  real nhpp_log(vector t, real beta, real theta, real tau){
    vector[num_elements(t)] loglik_part;
    real loglikelihood;
    for (i in 1:num_elements(t)){
      loglik_part[i] = log(beta) - beta*log(theta) + (beta - 1)*log(t[i]);
    }
    loglikelihood = sum(loglik_part) - (tau/theta)^beta;
    return loglikelihood;
  }
}
data {
  int<lower=0> N; //total # of failures
  int<lower=0> K; //number of predictors
  int<lower=0> S; //total # of shifts
  int<lower=0> D[S];//driver index, this must be an array
  vector<lower=0>[S] tau;//truncated time
  vector<lower=0>[N] event_time; //failure time
  int n_group[S]; //group sizes
  matrix[S, K] X_predictors;//predictor variable matrix
}
transformed data{
  matrix[S, K] X_centered;
  vector[S] X_means;
  for(m in 1:K){
    X_means[, m] = mean(X_predictors[, m]);
  }
}
parameters{
  real<lower=0> beta;
  vector[S] R0; // random intercept
  vector[3] R1_K; // fixed parameters
  real mu0; // hyperparameter
  real<lower=0> sigma0;// hyperparameter
}
transformed parameters{
  vector<lower=0>[S] theta;
  for (s0 in 1:S){
    theta[s0] = exp(R0[ D[s0] ] + X_predictors[s0,]*R1_K);
  }
}
model{
  int position;
  position = 1;
  for (s1 in 1:S){
    if(n_group[s1] == 0) continue;
    segment(event_time, position, n_group[s1]) ~ nhpp(beta, theta[s1], tau[s1]);
    position = position + n_group[s1];
  }
  beta ~ gamma(1, 1);
  R0 ~ normal(mu0, sigma0);
  R1_K  ~ normal(0, 10);
  mu0 ~ normal(0, 10);
  sigma0 ~ gamma(1, 1);
  theta ~ gamma(1, 0.01);
}
generated quantities{

}