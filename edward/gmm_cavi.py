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

N = 500  # number of data points
D = 2  # dimensionality of data

x_train = build_toy_dataset(N)

# Model
K = 2
pi = Dirichlet(tf.ones(K))
mu = Normal(tf.zeros(D), tf.ones(D), sample_shape=K)
sigmasq = InverseGamma(tf.ones(D), tf.ones(D), sample_shape=K)
x = ParamMixture(
    pi,
    {'loc': mu, 'scale_diag': tf.sqrt(sigmasq)},
    MultivariateNormalDiag,
    sample_shape=N
)
z = x.cat

z_t = ed.complete_conditional(z)
sess.run(z_t.prob([[0, 1], [2, 1]]))

sess.run(z.prob([0, 1]))
sess.run(z_t.logits)
p0 = z_t.prob(tf.ones([500, 1]))
