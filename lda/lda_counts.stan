data {
  int<lower=1> K; // num topics
  int<lower=1> V; // num words
  int<lower=0> D; // num docs
  int<lower=0> n[D, V]; // word counts for each doc

  // hyperparameters
  vector<lower=0>[K] alpha;
}

parameters {
  simplex[K] theta[D]; // topic mixtures
  simplex[V] beta[K]; // word dist for k^th topic
}

model {
  for (d in 1:D) {
    theta[d] ~ dirichlet(alpha);
  }

  for (d in 1:D) {
    vector[V] gamma;
    gamma = beta[1] * theta[d, 1];

    for (k in 2:K) {
      gamma = gamma + beta[k] * theta[d, k];
    }
    n[d] ~ multinomial(gamma);
  }

}
