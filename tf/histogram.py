"""

Playing around with histogram code from
https://www.tensorflow.org/get_started/tensorboard_histograms

date: 12/02/2017
"""

import tensorflow as tf

k = tf.placeholder(tf.float32)

# Make a normal distribution, with a shifting mean
mean_moving_normal = tf.random_normal(shape=[1000], mean=(5*k), stddev=1) ## this is a tensor, k will eventually be filled in with 1, ..., N
# Record that distribution into a histogram summary
# tf.summary is the general way of adding informatino to a tensorboard
tf.summary.histogram("normal/moving_mean", mean_moving_normal)

# Setup a session and summary writer
sess = tf.Session()
writer = tf.summary.FileWriter("/tmp/histogram_example")

summaries = tf.summary.merge_all()

# Setup a loop and write the summaries to disk
N = 400
for step in range(N):
  k_val = step/float(N)
  summ = sess.run(summaries, feed_dict={k: k_val})
  writer.add_summary(summ, global_step=step)
