"""

Experiment with Coordinate-Ascent Variational Inference (CAVI) for standard
Gaussian mixture models.

Combines ideas from
 - http://edwardlib.org/tutorials/unsupervised
 - http://edwardlib.org/api/ed/complete_conditional

author: sankaran.kris@gmail.com
date: 01/24/2018
"""

import numpy as np
import tensorflow as tf
import edward as ed
from edward.models import Dirichlet, MultivariateNormalDiag, Normal, \
  ParamMixture, PointMass

sess = tf.InteractiveSession()

# Data
def build_toy_dataset(N):
  pi = np.array([0.4, 0.6])
  mus = [[1, 1], [-1, -1]]
  stds = [[0.1, 0.1], [0.1, 0.1]]
  x = np.zeros((N, 2), dtype=np.float32)
  for n in range(N):
    k = np.argmax(np.random.multinomial(1, pi))
    x[n, :] = np.random.multivariate_normal(mus[k], np.diag(stds[k]))

  return x

N = 40  # number of data points
D = 2  # dimensionality of data

x_train = build_toy_dataset(N)

# Model
K = 3
pi = Dirichlet(tf.ones(K))
mu = Normal(tf.zeros(D), tf.ones(D), sample_shape=K)
sigmasq = InverseGamma(tf.ones(D), tf.ones(D), sample_shape=K)
x = []
z = []
for i in range(N):
  z.append(Categorical(probs=pi))
  x.append(Normal(mu[z[i].sample()], 1.0))

# Variational distributions
qz = []
for i in range(N):
  qz[i] = Categorical(tf.Variable(tf.random_dirichlet(K)))

qmu = Normal(
  loc=tf.Variable(tf.random_normal([D, K])),
  scale=tf.nn.softplus(tf.Variable(tf.random_normal([D, K])))
)

def local_updates(x_i, qmu):
  K = qmu.get_shape().as_list()[1]
  probs = np.zeros(K)

  for k in range(K):
    np.dot(qmu.mean()[:, k], x_i) - 0.5 * (qmu.variance()[:, K] + qmu.mean()[:, K]) ^ 2

# coordinate ascent
for iter in range(max_iter):

  for i in range(N):
    # update local variational params

  for k in range(K):
    # update global variational params

  # compute elbo
