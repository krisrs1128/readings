"""

Just following along the API description,
http://edwardlib.org/api/model

date: 01/01/2018
"""
import tensorflow as tf
from edward.models import Normal, Exponential, Dirichlet, MultivariateNormalTriL
tf.InteractiveSession()

# just defining some (vector and matrix) random variables
Normal(loc=tf.zeros(5), scale=tf.ones(5))
Exponential(rate=tf.ones([2, 3]))

# variable dimensions are the last coordinate
K = 10
Dirichlet(concentration=tf.constant([1.0] * K))
x = MultivariateNormalTriL(loc=tf.zeros([5, K]), scale_tril=tf.ones([5, K, K]))

x.eval()
x.log_prob(x).eval()
