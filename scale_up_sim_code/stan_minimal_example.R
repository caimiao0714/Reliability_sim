library(rstan)

stancode = 'data {
  int<lower=0> J;          // number of schools
  real y[J];               // estimated treatment effects
  real<lower=0> sigma[J];  // s.e. of effect estimates
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] eta;
}
transformed parameters {
  vector[J] theta;
  theta = mu + tau * eta;
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(y | theta, sigma);
}'

schools_data <- list(
  J = 8,
  y = c(28,  8, -3,  7, -1,  1, 18, 12),
  sigma = c(-15, 10, 16, 11,  9, 11, 10, 18)#Intentionally created a negative value here
)

for (i in 1:3) {
  tryCatch({fit1 <- stan(model_code  = stancode, data = schools_data,
                         chains = 1, iter = 1000, refresh = 0)}, error=function(e){})
}