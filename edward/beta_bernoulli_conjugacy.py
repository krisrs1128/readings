"""

Beta bernoulli example. Uses complete_conditionals methods. Following along
https://github.com/blei-lab/edward/blob/master/examples/beta_bernoulli_conjugate.py

date: 01/12/2018
"""

import edward as ed
import tensorflow as tf
import numpy as np

from edward.models import Bernoulli, Beta

sess = ed.get_session()

# Data
x_train = np.array([0, 1, 0, 0, 0, 0, 1, 1, 1, 0])

# Model
p = Beta(1., 1.)
x = Bernoulli(probs=p, sample_shape=10)

p_cond = ed.complete_conditional(p)

tf.global_variables_initializer().run()
alpha_bar = p_cond.parameters["concentration0"]
beta_bar = p_cond.parameters["concentration1"]
sess.run(alpha_bar, {x: x_train})
sess.run(beta_bar, {x: x_train})
