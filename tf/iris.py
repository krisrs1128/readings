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
