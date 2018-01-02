"""

Following along the inference API section,
http://edwardlib.org/api/inference-compositionality

date: 01/01/2017
"""
import edward as ed
from edward.models import Categorical, PointMass, Normal
import tensorflow as tf
sess = tf.Session()

# message passing example
N1 = 1000  # number of data points in first data set
N2 = 2000  # number of data points in second data set
D = 2  # data dimension
K = 5  # number of clusters

# MODEL
beta = Normal(loc=tf.zeros([K, D]), scale=tf.ones([K, D]))
z1 = Categorical(logits=tf.zeros([N1, K]))
z2 = Categorical(logits=tf.zeros([N2, K]))
x1 = Normal(loc=tf.gather(beta, z1), scale=tf.ones([N1, D]))
x2 = Normal(loc=tf.gather(beta, z2), scale=tf.ones([N2, D]))

# INFERENCE
qbeta = Normal(loc=tf.Variable(tf.zeros([K, D])),
               scale=tf.nn.softplus(tf.Variable(tf.zeros([K, D]))))
qz1 = Categorical(logits=tf.Variable(tf.zeros([N1, K])))
qz2 = Categorical(logits=tf.Variable(tf.zeros([N2, K])))

inference_z1 = ed.KLpq({beta: qbeta, z1: qz1}, {x1: x1_train})
inference_z2 = ed.KLpq({beta: qbeta, z2: qz2}, {x2: x2_train})

for _ in range(10000):
  inference_z1.update()
  inference_z2.update()
