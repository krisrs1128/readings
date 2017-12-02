"""

Following the tutorial at https://www.tensorflow.org/get_started/mnist/pros

author: sankaran.kris@gmail.com
date: 12/01/2017
"""

import tensorflow as tf
from tensorflow.examples.tutorials.mnist import input_data

###############################################################################
## functions used below
###############################################################################
def weight_variable(shape):
    """
    Weight initialization function
    """
    initial = tf.truncated_normal(shape, stddev=0.1)
    return tf.Variable(initial)


def bias_variable(shape):
    """Bias initialization function
    """
    initial = tf.constant(0.1, shape=shape)
    return tf.Variable(initial)


def conv2d(x, W):
    """Convolution wrapper
    """
    return tt.nn.conv2d(
        x,
        W,
        strides=[1, 1, 1, 1],
        padding="SAME"
    )


def max_pool_2x2(x):
    """Max pooling wrapper
    """
    return tf.nn.max_pool(
        x,
        ksize=[1, 2, 2, 1],
        strides=[1, 2, 2, 1],
        padding="SAME"
    )


###############################################################################
##  softmax regression
###############################################################################
mnist = input_data.read_data_sets('MNIST_data', one_hot=True)

x = tf.placeholder(
    tf.float32,
    shape = [None, 784] ## 28 x 28
)
y_ = tf.placeholder( ## presumably probabilities for 10 classes
    tf.float32,
    shape = [None, 10]
)

W = tf.Variable(tf.zeros([784, 10]))
b = tf.Variable(tf.zeros([10]))
y = tf.matmul(x, W) + b

cross_entropy = tf.reduce_mean(
    tf.nn.softmax_cross_entropy_with_logits(
        labels = y_,
        logits = y
    )
)

## initialize variables with random values
sess = tf.InteractiveSession()
sess.run(tf.global_variables_initializer())

## define optimization routine
train_step = tf.train.GradientDescentOptimizer(0.5).minimize(cross_entropy)
for _ in range(1000):
    batch = mnist.train.next_batch(100) ## batch size
    train_step.run(
        feed_dict = {x: batch[0], y_: batch[1]}
    )

## evaluate on test data
correct_pred = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
accuracy = tf.reduce_mean(
    tf.cast(correct_pred, tf.float32)
)
print(accuracy)
print(accuracy.eval(
    feed_dict={
        x: mnist.test.images,
        y_: mnist.test.labels
    }
))

###############################################################################
## Deeper net
###############################################################################
