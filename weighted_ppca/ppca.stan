/*
 * An example of standard PPCA
 *
 * author: sankaran.kris@gmail.com
 * date: 02/18/2017
 */


data {
  int < lower=0> N; // number of data points i n dataset
  int < lower=0> p; // dimension
  int < lower=0> k; // dimension of latent space
  vector [p] x[N]; // data
}
parameters {
  matrix[k, N] z;
  matrix[p, k] w;
  real<lower=0> sigma;
}
model {
  // priors
  to_vector (z) ~ normal(0,1);
  for (j in 1:p){
    w[j] ~ normal(0, 1);
  }

  // likelihood
  for (i in 1:N){
    x[i] ~ normal(w * col(z, i), sigma);
  }
}
