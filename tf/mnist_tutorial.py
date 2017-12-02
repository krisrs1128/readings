"""

Following the tutorial at https://www.tensorflow.org/get_started/mnist/pros

author: sankaran.kris@gmail.com
date: 12/01/2017
"""

from tensorflow.examples.tutorials.mnist import input_data
mnist = input_data.read_data_sets('MNIST_data', one_hot=True)

x = tf.placeholder(
    tf.float32,
    shape = [None, 784] ## 28 x 28
)
y = tf.placeholder( ## presumably probabilities for 10 classes
    tf.float32,
    shape = [None, 10]
)

W = tf.Variable(tf.zeros([784, 10]))
b = tf.Variable(tf.zeros([10]))
y = tf.matmult(x, W) + b

cross_entropy = tf.reduce_mean(
    tf.nn.softmax_cross_entropy_with_logits(
        labels = y_,
        logits = y
    )
)

## define optimization routine
train_step = tf.train.GradientDescentOptimizer(0.5).minimize(cross_entropy)
for _ in range(1000):
    batch = mnist.train.next_batch(100) ## batch size
    train_step.run(
        feed_dict = {x: batch[0], y_: batch[1]}
    )


## initialize variables with random values
sess.run(tf.global_variables_initializer())
