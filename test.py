import numpy as np
import inspect as insp
import math
import scipy.linalg as sp
import threeD_rotation as td
import Lorentz_Transformation as lt
import LorentzGenerators as lg

# for i in range(0,100):
#     a = np.random.rand(4,)
#     b = lt.Lorentz4vector(components=a,name='a',mass=1)
#     b.lorentz_rot_toz()
#     c = b.get_4_vector()
#     print(c[2]+c[1])

#test passed

vector_1 = [85,84,12,4]
vector_2 = [2,1,1,1]
v1 = lt.Lorentz4vector(components=vector_1)
v2 = lt.Lorentz4vector(components=vector_2)
for i in range(100):
    vector = np.random.rand(4)
    vector[0] += 1
    v3 = lt.Lorentz4vector(components=vector)
    if v3.get_on_shell_ness():
        v3.lorentz_rot_toz()
        v3.zrapidity()
        v3.z_boost_to_rest_frame()
        print(np.sum(v3.get_4_vector()[1:]))
# print(lt.four_vector_contraction_component(vector_1,vector_2))
# print(lt.general_matrix_form(parameter=[np.pi/2,0,0,0,0,0]))