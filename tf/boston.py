"""

Example from Boston Housing ata set. Used to explain input_funs. This is
basically a way of preparing the training and test data for subsequent
modeling.

https://www.tensorflow.org/get_started/input_fn

date: 12/02/2017
"""

import itertools
import pandas as pd
import tensorflow as tf
tf.logging.set_verbosity(tf.logging.INFO)

###############################################################################
## functions
###############################################################################
def get_input_fn(data_set, num_epochs=None, shuffle=True):
    """Input Fun. Functional

    Let's you use the same code for train, validation, and test
    """
    return tf.estimator.inputs.pandas_input_fn(
        x=pd.DataFrame({k: data_set[k].values for k in FEATURES}),
        y = pd.Series(data_set[LABEL].values),
        num_epochs=num_epochs,
        shuffle=shuffle)

###############################################################################
## read in data
###############################################################################
COLUMNS = ["crim", "zn", "indus", "nox", "rm", "age",
           "dis", "tax", "ptratio", "medv"]
FEATURES = ["crim", "zn", "indus", "nox", "rm",
            "age", "dis", "tax", "ptratio"]
LABEL = "medv"

training_set = pd.read_csv("boston_train.csv", skipinitialspace=True,
                           skiprows=1, names=COLUMNS)
test_set = pd.read_csv("boston_test.csv", skipinitialspace=True,
                       skiprows=1, names=COLUMNS)
prediction_set = pd.read_csv("boston_predict.csv", skipinitialspace=True,
                             skiprows=1, names=COLUMNS)

# defines feature columns (more useful in case we have mixed types)
feature_cols = [tf.feature_column.numeric_column(k) for k in FEATURES]
regressor = tf.estimator.DNNRegressor(
    feature_columns=feature_cols,
    hidden_units=[10, 10],
    model_dir="/tmp/boston_model"
)

regressor.train(input_fn=get_input_fn(training_set), steps=5000)

###############################################################################
## evaluation
###############################################################################
ev = regressor.evaluate(
    input_fn=get_input_fn(test_set, num_epochs=1, shuffle=False
    ))
loss_score = ev["loss"]
print("Loss: {0:f}".format(loss_score))
