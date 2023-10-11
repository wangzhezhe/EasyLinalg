import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# generate some data from a multivariate gaussin distribution 
# which have three variables
# mu=np.array([1,10,20])
# sigma=np.matrix([[4,10,0],[10,25,0],[0,0,100]])
# data=np.random.multivariate_normal(mu,sigma,1000)
# values = data.T

# # size of the values are 3*1000
# # three attribute, each has 1000 samples
# print(values.shape)

# kde = stats.gaussian_kde(values)

# # density = kde(values)
# # fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
# # x, y, z = values
# # ax.scatter(x, y, z, c=density)
# # plt.savefig('kde_test.png')

# # test value [0,1,2]
# density = kde([[0,1,2,3],[3,4,5,1],[1,3,2,4]])

# print("density",density)

# test kde2 

test_data = np.matrix([[0,0.15,-0.1],[2.0,2.3,2.5]])

print(test_data.shape)

kde = stats.gaussian_kde(test_data)
print("evaluate pdf", kde.evaluate([0.1,2.5]))
