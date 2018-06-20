import sympy as sym
import scipy.linalg as sp
import threeD_rotation as td
import numpy as np
import Lorentz_Transformation as lt
import LorentzGenerators as lg
import MonteCarlo as mc
import time as time

test_decay0 = lt.Lorentz4vector(name='test',components=[5,0,0,3],mass=4)
test_decay1 = lt.Lorentz4vector(name='test',components=[5,0,0,3],mass=4)
test_decay3 = lt.Lorentz4vector(components=[5,0,0,3],mass=4)


angles = np.random.rand(2,)
print(np.cos(angles[0]*np.pi))
mat1 = lt.general_matrix_form(parameter=[0,np.pi*angles[0],0,0,0,0])
test_decay0.lorentz_transformation_off_matrix(mat1,mode=0)
# print(test_decay0)
print(lt.cos_theta(test_decay0,test_decay1))
mat2 = lt.general_matrix_form(parameter=[0,0,2*np.pi*angles[1],0,0,0])
test_decay0.lorentz_transformation_off_matrix(mat2)
print(lt.cos_theta(test_decay0,test_decay1))
mat = np.matmul(mat2,mat1)
test_decay3.lorentz_transformation_off_matrix(mat)
print(lt.cos_theta(test_decay3,test_decay1))
print(test_decay0,test_decay3)

# Surprisingly, the originial naive algorithm seems the right

# mat2 = lt.general_matrix_form(parameter=[0,0,2*np.pi*angles[1],0,0,0])
# mat = np.matmul(mat1,mat2)
mat = np.diag((1,np.sin(angles[0])*np.cos(angles[1]),np.sin(angles[0])*np.sin(angles[1]),np.cos(angles[0])))
# print(mat)
# Sadly, this dumb way works better than the sophisticated generator way.

test_decay0.lorentz_transformation_off_matrix(mat)
# print(lt.cos_theta(test_decay0,test_decay1))
# print(angles)


# print("Running ti me : %f CPU seconds" % (end_CPU - start_CPU))
