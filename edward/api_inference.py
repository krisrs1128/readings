"""

Following along the inference API section,
http://edwardlib.org/api/inference

date: 01/01/2017
"""
import edward as ed
from edward.models import Beta, Bernoulli
import tensorflow as tf
sess = tf.Session()

# Beta-Bernoulli model
pi = Beta(1.0, 1.0)
x = Bernoulli(probs=pi, sample_shape=10)

# Beta posterior; it conditions on the sample tensor associated to x
pi_cond = ed.complete_conditional(pi)

# Generate samples from p(pi | x = NumPy array)
sess.run(pi_cond, {x: np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 1])})
