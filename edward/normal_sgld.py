"""

Correlated normal posterior sampling, following
https://github.com/blei-lab/edward/blob/master/examples/normal_sgld.py

date: 01/12/2018
"""

import edward as ed
import tensorflow as tf

from edward.models import MultivariateNormalTriL, Empirical

# multivariate normal model
cov = [[ 0.36,  0.12,  0.06],
       [ 0.12,  0.29, -0.13],
       [ 0.06, -0.13,  0.26]]

z = MultivariateNormalTriL(
    [1., 2, 3],
    tf.cholesky(cov)
)

x = [z.eval() for i in range(100)]
x = np.array(x)

qz = Empirical(params=tf.Variable(tf.random_normal([2000, 3])))
inference = ed.SGLD({z: qz})
inference.run(step_size=0.5)

sess = ed.get_session()
qz.value().eval()
