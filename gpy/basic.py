"""

Following along tutorial at
http://nbviewer.jupyter.org/github/SheffieldML/notebook/blob/master/GPy/models_basic.ipynb

author: sankaran.kris@gmail.com
date: //2018
"""

#import necessary modules, set up the plotting
import numpy as np
from matplotlib import pyplot as plt
import GPy

m = GPy.examples.regression.sparse_GP_regression_1D(plot=False, optimize=False)

print(m)
print(m.rbf)
print(m.inducing_inputs)

m.plot()
plt.show()
plt.close()

m.rbf.lengthscale = 0.1
m.plot()
plt.show()
plt.close()

m.optimize()
m.plot()
plt.show()
plt.close()

print(m.gradient)

m.plot(plot_density=True)
plt.show()

# some simulation
# http://gpss.cc/gpss13/assets/lab1.pdf
k = GPy.kern.RBF(1, 1., 0.2)
n = 100
x = np.random.random([n, 1])
mu = np.zeros(n)
C = k.K(x, x) + 0.01 * np.eye(n)
y = np.random.multivariate_normal(np.zeros(n), C)
plt.scatter(x[:], y.transpose())
plt.show()
