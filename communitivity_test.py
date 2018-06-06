import Lorentz_Transformation as lt
import LorentzGenerators as lg
import numpy as np

p1 = lt.Lorentz4vector(components=[5,0,0,5],mass=0)
p2 = lt.Lorentz4vector(components=[5,0,0,-5],mass=0)

angles = np.random.rand(4,)

mat1 = lt.general_matrix_form(parameter=[0,0,angles[0],0,0,0])
mat2 = lt.general_matrix_form(parameter=[0,angles[1],0,0,0,0])

mat3 = lt.general_matrix_form(parameter=[0,0,angles[2],0,0,0])
mat4 = lt.general_matrix_form(parameter=[0,angles[3],0,0,0,0])

print(np.matmul(mat1,mat3)-np.matmul(mat3,mat1))


c = lt.general_matrix_form(parameter=[0,angles[1],angles[0],0,0,0])
a = np.matmul(mat1,mat2)
b = np.matmul(mat2,mat1)
# print(a)
# print(b)
# print(np.matmul(a,p1.components))
# print(np.matmul(b,p1.components))
# print(np.matmul(c,p1.components))

####################################################
# Conclusions:
# Generaly rotations do NOT commute with each other.
# Rotation with respect to the same axis commutes with each other.