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

K4_2=np.array(
 [
    [23.1819286,21.6486893,21.1437302,26.6601410],
    [21.6486893,38.5402031,20.2783680,6.3069100],
    [21.1437302,20.2783680,20.0003567,23.6239624],
    [26.6601410,6.3069100,23.6239624,50.3645554]
 ]
)

K4_3=np.array(
    [
        [1.0, 0.9, 0.0, 0.0],
        [0.9, 1.0, 0.9, 0.0],
        [0.0, 0.9, 1.0, 0.9],
        [0.0, 0.0, 0.9, 1.0]
    ]
)


# w, v = LA.eig(K4_1)

# print(w)
# print(v)

# w, v = LA.eig(K4_2)

# print(w)
# print(v)


w, v = LA.eig(K4_3)

print(w)
print(v)