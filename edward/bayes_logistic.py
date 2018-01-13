"""

Following along the example from
https://github.com/blei-lab/edward/blob/master/examples/bayesian_logistic_regression.py

date: 01/12/2018
"""

import edward as ed
import numpy as np
import tensorflow as tf

from edward.models import Bernoulli, Normal, Empirical

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def build_toy_data(N, D=1, sigma=0.1):
    X = np.random.normal(size=(N, D))
    beta = np.random.normal(size=(D, 1))
    eta = sigmoid(np.dot(X, beta) + sigma * np.ones((N, 1)))
    y = 1. * (np.random.randn(N, 1) > eta)

    return X, y


N = 100
D = 1
X_train, y_train = build_toy_data(N, D)

# Model
X = tf.placeholder("float32", [N, D])
sigma_beta = 10.
beta = Normal(tf.zeros((D, 1)), sigma_beta)
y = Bernoulli(tf.matmul(X, beta))

# (Variational) Inference
qloc = tf.Variable(tf.random_normal((D, 1)))
qscale = tf.Variable(tf.random_gamma((D, 1), 1))
qbeta = Normal(loc=qloc, scale=qscale)

inference = ed.KLqp({beta: qbeta}, data = {X: X_train, y: y_train})
inference.initialize(n_iter=500)
