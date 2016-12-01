/*
 * Dynamic Topic Model
 *
 * This models documents using a sequence of cluster and topic
 * parameters that evolve over time. See Blei & Lafferty 2006.
 */

data {
  int<lower=0> N; // number of samples
  int<lower=0> V; // number of words
  int<lower=0> T; // number of unique times
  int<lower=0> K; // number of topics
  real<lower=0> sigma; // rate of topic evolution on simplex
  real<lower=0> delta; // rate of membership evolution on simplex

  real times[T]; // unique times
  int<lower=0> times_mapping[N]; // times associated to each sample
  int<lower=0> X[N, V]; // word counts
}

parameters {
  vector[V] mu[T, K]; // logitted vocabulary probabilities
  vector[K] alpha[T]; // logitted cluster probabilities
}

transformed parameters {
  vector[V] beta[T, K];
  vector[K] theta[T];
  for (i in 1:T) {
    for (k in 1:K) {
      beta[i, k] = softmax(mu[i, k]);
    }
    theta[i] = softmax(alpha[i]);
  }
}

model {
  for (i in 1:(T - 1)) {
    for (k in 1:K) {
      mu[i + 1, k] ~ normal(mu[i, k], sqrt(times[i + 1] - times[i]) * sigma);
    }
    alpha[i + 1] ~ normal(alpha[i], sqrt(times[i + 1] - times[i]) * delta);
  }

  for (i in 1:N) {
    vector[V] gamma;
    gamma = beta[times_mapping[i], 1] * theta[times_mapping[i]][1];
    for (k in 2:K) {
      gamma = gamma + beta[times_mapping[i], k] * theta[times_mapping[i]][k];
    }

    X[i] ~ multinomial(gamma);
  }
}
