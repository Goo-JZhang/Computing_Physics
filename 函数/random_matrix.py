import numpy as np

rd = np.random.RandomState(0)
A = np.zeros((20,6,6))
for i in range(20):
    A[i] = rd.uniform(-1,1,(6,6))
np.savetxt("random_matrix_py.txt",A[0])