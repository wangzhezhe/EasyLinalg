import numpy as np 
import matplotlib.pyplot as plt
from numpy import linalg as LA
from sklearn.datasets import make_spd_matrix


K4_1=np.array(
    [
        [0.20, 0.60, 0.40, 0.80],
        [0.60, 1.80, 1.20, 2.40],
        [0.40, 1.20, 0.80, 1.60],
        [0.80, 2.40, 1.60, 3.20]
    ]
)


w, v = LA.eig(K4_1)

print(w)
print(v)