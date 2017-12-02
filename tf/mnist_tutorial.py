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

    see https://www.tensorflow.org/api_docs/python/tf/nn/conv2d
    """
    return tf.nn.conv2d(
        x,
        W, ## this is the filter applied to each patch
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
x_image = tf.reshape(x, [-1, 28, 28, 1])

W_conv1 = weight_variable([5, 5, 1, 32]) ## 1 input, 32 output channels
b_conv1 = bias_variable([32])
h_conv1 = tf.nn.relu(conv2d(x_image, W_conv1) + b_conv1)
h_pool1 = max_pool_2x2(h_conv1)

## second layer
W_conv2 = weight_variable([5, 5, 32, 64])
b_conv2 = bias_variable([64])
h_conv2 = tf.nn.relu(
    conv2d(h_pool1, W_conv2) + b_conv2
)
h_pool2 = max_pool_2x2(h_conv2)

## last layer
W_fc1 = weight_variable([7 * 7 * 64, 1024])
b_fc1 = bias_variable([1024])
h_pool2_flat = tf.reshape(h_pool2, [-1, 7 * 7 * 64])
h_fc1 = tf.nn.relu(tf.matmul(h_pool2_flat, W_fc1) + b_fc1)

## dropout
keep_prob = tf.placeholder(tf.float32)
h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

W_fc2 = weight_variable([1024, 10])
b_fc2 = bias_variable([10])
y_conv = tf.matmul(h_fc1_drop, W_fc2) + b_fc2

## specify training for this model
cross_entropy = tf.reduce_mean(
    tf.nn.softmax_cross_entropy_with_logits(
        labels=y_,
        logits=y_conv
    )
)
train_step = tf.train.AdamOptimizer(1e-4).minimize(cross_entropy)
correct_pred = tf.equal(
    tf.argmax(y_conv, 1),
    tf.argmax(y_, 1)
)
accuracy = tf.reduce_mean(
    tf.cast(correct_pred, tf.float32)
)

## train the model
with tf.Session() as sess:
    sess.run(tf.global_variables_initializer())
    for i in range(20000):
        batch = mnist.train.next_batch(50)
        if i % 20 == 0:
            train_acc = accuracy.eval(
                feed_dict={
                    x: batch[0],
                    y_: batch[1],
                    keep_prob: 1.0
                }
            )
            print("i %d accuracy %g" % (i, train_acc))
        train_step.run(
            feed_dict={
                x: batch[0],
                y_: batch[1],
                keep_prob: 0.5
            })


## test accuracy
acc = accuracy.eval(
    feed_dict={
        x: mnist.test.images,
        y_: mnist.test.labels,
        keep_prob: 1.0
    }
)
