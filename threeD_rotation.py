import numpy as np
import scipy.linalg as sp

so3g = np.zeros((3,3,3))
so3g[0,2,1]=so3g[1,0,2]=so3g[2,1,0]=-1
so3g[0,1,2]=so3g[1,2,0]=so3g[2,0,1]=1
# print(so3g)

#This does not work as expected. Abandoned!
def rotation(vector,angle):
    summ = 0
    for i in range(3):
        summ -= angle[i]*so3g[i]
    # print(summ)
    rotation_matrix = sp.expm(summ)
    # print(rotation_matrix)
    vector_1 = np.matmul(rotation_matrix,vector)
    return vector_1

#Find the correct rotation angle to rotate a certain 3D vector to z direction, or the opposite
#To specific a spatical angle one needs only two variables: theta and phi
#We first rotate the vector into the xz plane(phi rotation), then into z direction(theta rotation)

def angles(vector):
    norm = np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
    phi = np.arccos(vector[0]/np.sqrt(vector[0]**2+vector[1]**2))
    theta = np.arccos(vector[2]/norm)
    # print(phi,theta)
    return [phi,theta]

def rotation_matrix(phi,theta):
    mat_rotphi = sp.expm(phi*so3g[2])
    mat_rottheta = sp.expm(theta*so3g[1])
    matrix = np.matmul(mat_rottheta,mat_rotphi)
    return matrix

# print(angles([7,24,25]))
# print(np.matmul(rotation_matrix(*angles([5,12,13])),[5,12,13]))
#     try :
#         cos_x = vector[0]/np.sqrt(vector[1]**2+vector[0]**2)
#         angle_x = np.arccos(cos_x)
#     except ValueError:
#         angle_x = 0
#     try:
#         cos_y = vector[0]/np.sqrt(vector[1]**2+vector[0]**2)
#         angle_y = np.arccos(vector[1]/np.sqrt(vector[2]**2+vector[1]**2))
#     except ValueError:
#         angle_y = 0
#     try:
#         angle_z = np.arccos(vector[2]/np.sqrt(vector[2]**2+vector[0]**2))
#     except ValueError:
#         angle_z = 0
#     return [angle_x,angle_y,angle_z]
# vector = np.random.rand(3,)
# alpha = angles(vector)
# print(rotation(vector,angle=[alpha[0],alpha[1],0]))


#Key function of this module. 3D rotate an arbitrary vectro into z direction
def rotate_to_z_axis(vector):
    return np.matmul(rotation_matrix(*angles(vector)),vector)

def norm(vector):
    return np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
# print(rotate_to_z_axis([3,4,5]))
for i in range(100):
    vector = np.random.rand(3)
    if norm(vector) - norm(rotate_to_z_axis(vector)) > 1e-8:
        print("Alert!")
    count = 0
    for j in rotate_to_z_axis(vector):
        if j > 1e-12:
            count += 1
    if count > 1:
        print("Alert!")

#Test Passed!
