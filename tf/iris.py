"""

Following along the iris example for learning the tf.estimator API

https://www.tensorflow.org/get_started/estimator

date: 12/02/2017
"""

import os
from six.moves.urllib.request import urlopen

import tensorflow as tf
import numpy as np

###############################################################################
## Get the data
###############################################################################
IRIS_TRAINING = "iris_training.csv"
IRIS_TRAINING_URL = "http://download.tensorflow.org/data/iris_training.csv"

IRIS_TEST = "iris_test.csv"
IRIS_TEST_URL = "http://download.tensorflow.org/data/iris_test.csv"

if not os.path.exists(IRIS_TRAINING):
  raw = urlopen(IRIS_TRAINING_URL).read()
  with open(IRIS_TRAINING,'wb') as f:
    f.write(raw)

if not os.path.exists(IRIS_TEST):
  raw = urlopen(IRIS_TEST_URL).read()
  with open(IRIS_TEST,'wb') as f:
    f.write(raw)

# Load datasets.
training_set = tf.contrib.learn.datasets.base.load_csv_with_header(
    filename=IRIS_TRAINING,
    target_dtype=np.int,
    features_dtype=np.float32)
test_set = tf.contrib.learn.datasets.base.load_csv_with_header(
    filename=IRIS_TEST,
    target_dtype=np.int,
    features_dtype=np.float32)

###############################################################################
##  specify the model and train
###############################################################################
# Specify that all features have real-value data
feature_columns = [tf.feature_column.numeric_column("x", shape=[4])]

# Build 3 layer DNN with 10, 20, 10 units respectively.
classifier = tf.estimator.DNNClassifier(
    feature_columns=feature_columns,
    hidden_units=[10, 20, 10],
    n_classes=3,
    model_dir="/tmp/iris_model"
)

# Define the training inputs
# before training, need to specify how to process input data
train_input_fn = tf.estimator.inputs.numpy_input_fn(
    x={"x": np.array(training_set.data)},
    y=np.array(training_set.target),
    num_epochs=None,
    shuffle=True)

# Train model.
classifier.train(input_fn=train_input_fn, steps=2000)

###############################################################################
## Evaluation
###############################################################################
# Define the test inputs
test_input_fn = tf.estimator.inputs.numpy_input_fn(
    x={"x": np.array(test_set.data)},
    y=np.array(test_set.target),
    num_epochs=1,
    shuffle=False)

# Evaluate accuracy.
accuracy_score = classifier.evaluate(input_fn=test_input_fn)["accuracy"]
