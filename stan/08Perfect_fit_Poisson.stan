data {
  int N;
  int K;
  int Y[N];
  int id[N];
  real x1[N];
  real x2[N];
}
parameters{
  real mu0;
  real sigma0;
  vector[K] b0;
  real b1;
  real b2;
}
model {
  // intermediate parameters;
  vector[N] tmp;
  for(i in 1:N){
    tmp[i] = b0[id[i]] + b1*x1[i] + b2*x2[i];
  }

  // increase log likelihood;
  target += poisson_lpmf(Y|tmp);

  //priors
  b0 ~ normal(mu0, sigma0);
  mu0 ~ normal(0, 5);
  sigma0 ~ gamma(1, 1);
  b1 ~ normal(0, 10);
  b2 ~ normal(0, 10);
}
