/*
 * Dynamic unigram model
 *
 * This models documents whose word proportions evolve according to a
 * softmax on an ordinary random walk. This is the one-topic version
 * of the kalman-filter based dynamic topic model.
 */

data {
  int<lower=0> N; // number of samples
  int<lower=0> V; // number of words
  int<lower=0> T; // number of unique times
  real<lower=0> sigma; // rate of evolution on simplex

  real times[T]; // unique times
  int<lower=0> times_mapping[N]; // times associated to each sample
  int<lower=0> X[N, V]; // word counts
}

parameters {
  vector[V] beta[T];
}

model {
  for (i in 1:(T - 1)) {
    beta[i + 1] ~ normal(beta[i], sqrt(times[i + 1] - times[i]) * sigma);
  }

  for (i in 1:N) {
    X[i] ~ multinomial(softmax(beta[times_mapping[i]]));
  }
}
