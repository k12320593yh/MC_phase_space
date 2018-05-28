import numpy as np
# import pandas as p
minkowski_metric = np.array([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]])
lorentz_generators = [0,0,0,0,0,0,0]
base = np.zeros((4,4))
rotation_1,rotation_2,rotation_3,boost_1,boost_2,boost_3 = np.zeros((4,4)),np.zeros((4,4)),np.zeros((4,4)),np.zeros((4,4)),np.zeros((4,4)),np.zeros((4,4))
rotation_1[2,3] = -1
rotation_2[3,1] = -1
rotation_3[1,2] = -1
rotation_1[3,2] = 1
rotation_2[1,3] = 1
rotation_3[2,1] = 1
boost_1[1,0] = boost_1[0,1] = 1
boost_2[0,2] = boost_2[2,0] = boost_3[0,3] = boost_3[3,0] = 1
identity = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
lorentz_generators[0] = identity
lorentz_generators[1] = rotation_1
lorentz_generators[2] = rotation_2
lorentz_generators[3] = rotation_3
lorentz_generators[4] = boost_1
lorentz_generators[5] = boost_2
lorentz_generators[6] = boost_3
lorentz_generators = np.array(lorentz_generators)
gamma_0 = np.array([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])
gamma_1 = np.array([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]])
gamma_2 = np.array([[0,0,0,-1j],[0,0,1j,0],[0,1j,0,0],[-1j,0,0,0]])
gamma_3 = np.array([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]])
gamma_5 = 1j*np.matmul(gamma_0,np.matmul(gamma_1,np.matmul(gamma_2,gamma_3)))
# print(gamma_5)
gamma_matrices = np.array([gamma_0,gamma_1,gamma_2,gamma_3])
# print(gamma_matrices[3])

# for i in range(4):
#     print(gamma_matrices[i])
# print(gamma_2.imag)
# # for i in range(7):
#     print(lorentz_generators[i])