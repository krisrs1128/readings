"""

Following the prob PCA discussion

Mainly interested in how you'd evaluate predictive likelihood.
http://edwardlib.org/tutorials/probabilistic-pca

date: 01/02/2018
"""

import numpy as np
import edward as ed
import tensorflow as tf
from edward.models import Normal

tf.InteractiveSession()

########################################
# Utils
########################################

def build_toy_dataset(N, D, K, sigma=1):
    x_train = np.zeros((D, N), dtype="float32")
    w = np.random.normal(0.0, 2.0, size=(D, K)).astype("f")
    z = np.random.normal(0.0, 1.0, size=(K, N))
    mean = np.dot(w, z)
    for d in range(D):
        for n in range(N):
            x_train[d, n] = np.random.normal(mean[d, n], sigma)

    return w, x_train


########################################
# Generate data
########################################
N = 5000  # number of data points
D = 2  # data dimensionality
K = 1  # latent dimensionality

w_true, x_train = build_toy_dataset(N, D, K)

########################################
# Setup model and inference
########################################
w = Normal(loc=tf.zeros([D, K]), scale=2.0 * tf.ones([D, K]))
z = Normal(loc=tf.zeros([N, K]), scale=tf.ones([N, K]))
x = Normal(loc=tf.matmul(w, z, transpose_b=True), scale=tf.ones([D, N]))

qw = Normal(loc=tf.Variable(tf.random_normal([D, K])),
            scale=tf.nn.softplus(tf.Variable(tf.random_normal([D, K]))))
qz = Normal(loc=tf.Variable(tf.random_normal([N, K])),
            scale=tf.nn.softplus(tf.Variable(tf.random_normal([N, K]))))

inference = ed.KLqp({w: qw, z: qz}, data={x: x_train})
inference.run(n_iter=500, n_print=100, n_samples=10)

print(qw.mean().eval())

x_post = ed.copy(x, {z: qz, w: qw})
ed.evaluate("log_likelihood", {x_post: x_train})

x_unrelated = np.random.normal(size=(2, 5000)).astype("f")
ed.evaluate("log_likelihood", {x_post: x_unrelated})
