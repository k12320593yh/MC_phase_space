import numpy as np
import inspect as insp
import math
import scipy.linalg as sp
import threeD_rotation as td
import Lorentz_Transformation as lt
import LorentzGenerators as lg
import MonteCarlo as mc

weight_all_zero = 6.332573977646112e-09
tw = 0.2223
mz = 91.2
mmu = 105.7
me = 0.5
mtau = 1776.8
mw = 80.3
e = np.sqrt(4*np.pi/137)
# print(e)
gz = 0.718
gw = 0.629
alpha = 1/137
gf = 1.11639e-5

# for i in range(0,100):
#     a = np.random.rand(4,)
#     b = lt.Lorentz4vector(components=a,name='a',mass=1)
#     b.lorentz_rot_toz()
#     c = b.get_4_vector()
#     print(c[2]+c[1])

#test passed

# vector_1 = [85,84,12,4]
# vector_2 = [2,1,1,1]
# v1 = lt.Lorentz4vector(components=vector_1)
# v2 = lt.Lorentz4vector(components=vector_2)
# for i in range(100):
#     vector = np.random.rand(4)
#     vector[0] += 1
#     v3 = lt.Lorentz4vector(components=vector)
#     if v3.get_on_shell_ness():
#         v3.lorentz_rot_toz()
#         v3.zrapidity()
#         v3.z_boost_to_rest_frame()
#         print(np.sum(v3.get_4_vector()[1:]))
# print(lt.four_vector_contraction_component(vector_1,vector_2))
# print(lt.general_matrix_form(parameter=[np.pi/2,0,0,0,0,0]))

def test_amplitude(vectors,masses):
    # print(vectors.shape)
    a = lt.fcc(vectors[0],vectors[1])
    return a

def test_amplitude_3():
    p = 0
    p = p+1


def eenunugamma(vectors, mz=91.2, p=None,masses=[0,0,0]):
    p = np.vstack((np.array([[0,0,0,0],[5,0,0,5],[5,0,0,-5]]),vectors))
    # print(p)
    coefficient = e**2*gz**4
    b = (16*lt.fcc(p[1],p[5])*lt.fcc(p[2],p[5])*(mz**2-2*(lt.fcc(p[3],p[4])))**2)
    a = (4*tw**2*lt.fcc(p[1],p[3])*lt.fcc(p[1],p[5])*lt.fcc(p[2],p[4])-
         8*tw**2*lt.fcc(p[1],p[5])*lt.fcc(p[2],p[3])*lt.fcc(p[2],p[4])+
         4*tw**2*lt.fcc(p[1],p[3])*lt.fcc(p[2],p[4])*lt.fcc(p[2],p[5])+
         (lt.fcc(p[1],p[5]))*(4*tw**2*lt.fcc(p[1],p[5])*(lt.fcc(p[2],p[3])+lt.fcc(p[3],p[5])) -lt.fcc(p[2],p[5])*(8*tw**2-4*tw+1)*(lt.fcc(p[1],p[3]))-
         4*tw**2*lt.fcc(p[2],p[3]))+
         lt.fcc(p[1],p[2])*(4*tw**2*lt.fcc(p[1],p[4])*(lt.fcc(p[3],p[5])-2*lt.fcc(p[2],p[3])))+
         (-1+2*tw**2)*lt.fcc(p[1],p[3])*(2*lt.fcc(p[2],p[4])-lt.fcc(p[4],p[5]))+
         4*tw**2*(lt.fcc(p[2],p[4])*lt.fcc(p[3],p[5])+lt.fcc(p[2],p[3])*lt.fcc(p[4],p[5]))-(4*tw-1)*lt.fcc(p[2],p[4])*lt.fcc(p[3],p[5])-
         4*tw*lt.fcc(p[1],p[3])*lt.fcc(p[1],p[5])*lt.fcc(p[2],p[4])+
         4*tw*(lt.fcc(p[1],p[5])*lt.fcc(p[2],p[3])*lt.fcc(p[2],p[4])-lt.fcc(p[1],p[3])*lt.fcc(p[2],p[4])*lt.fcc(p[2],p[5]))+
         lt.fcc(p[1],p[3])*lt.fcc(p[1],p[5])*lt.fcc(p[2],p[4])-lt.fcc(p[1],p[3])*lt.fcc(p[2],p[4])*lt.fcc(p[2],p[5])+
         lt.fcc(p[1],p[3])*lt.fcc(p[2],p[4])*lt.fcc(p[2],p[5])+4*tw**2*lt.fcc(p[1],p[3])*lt.fcc(p[1],p[5])*lt.fcc(p[4],p[5])-
         (4*tw+1)*lt.fcc(p[1],p[3])*lt.fcc(p[1],p[5])*lt.fcc(p[4],p[5])+4*tw**2*(lt.fcc(p[2],p[4])*lt.fcc(p[2],p[5])*lt.fcc(p[3],p[5])+
         lt.fcc(p[2],p[3])*lt.fcc(p[2],p[5])*lt.fcc(p[4],p[5]))-
         (4*tw-1)*lt.fcc(p[2],p[4])*lt.fcc(p[2],p[5])*lt.fcc(p[3],p[5]))
    # print(a)
    # print(b)
    return coefficient*a/b

def eenunugamma1(vectors,masses,s_sqrt=500):
    p1 = [s_sqrt/2,0,0,s_sqrt/2]
    p2 = [s_sqrt/2,0,0,(0 - 1)*s_sqrt/2]
    k1 = vectors[0]
    k2 = vectors[1]
    k3 = vectors[2]
    # print(p1,p2,k1,k2,k3)
    # print(k1+k2+k3)
    sw2 = 0.222246494
    cw2 = 1-sw2
    mz2 = 91.2**2
    a = (e**6*((-16*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (64*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (128*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (16*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (64*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (128*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (32*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (128*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (256*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (32*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (128*sw2*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (256*sw2**2*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (32*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (128*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (256*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (16*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (64*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (128*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (32*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (128*sw2*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (256*sw2**2*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (16*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (64*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (128*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + 2*(-4/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (16*sw2)/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2)*(lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (p1)) - lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)) - lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p1)) + lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)) + lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p1), (p2)) - lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p1), (p2))) - (16*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (64*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (128*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (16*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (64*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (128*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (32*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (128*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (256*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (32*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (128*sw2*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (256*sw2**2*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (32*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (128*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (256*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (32*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (128*sw2*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (256*sw2**2*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (16*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (64*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (128*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (16*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (64*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (128*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (32*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (32*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (64*lt.fcc((k1), (p2))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (512*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (64*lt.fcc((k1), (p1))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (512*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (32*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (32*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (32*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (128*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (256*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (32*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (128*sw2*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (256*sw2**2*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (64*lt.fcc((k1), (p2))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (256*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (512*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (p2))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (32*lt.fcc((k1), (k2))*lt.fcc((k3), (p2))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (k2))*lt.fcc((k3), (p2))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (k2))*lt.fcc((k3), (p2))*lt.fcc((p1), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (32*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (32*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (32*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (64*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (512*sw2**2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (32*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (64*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (512*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (32*lt.fcc((k1), (k2))*lt.fcc((k3), (k3))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (128*sw2*lt.fcc((k1), (k2))*lt.fcc((k3), (k3))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (256*sw2**2*lt.fcc((k1), (k2))*lt.fcc((k3), (k3))*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (32*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (128*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (256*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (32*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (128*sw2*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (256*sw2**2*lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (64*lt.fcc((k1), (p1))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (256*sw2*lt.fcc((k1), (p1))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (512*sw2**2*lt.fcc((k1), (p1))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (32*lt.fcc((k1), (k2))*lt.fcc((k3), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (128*sw2*lt.fcc((k1), (k2))*lt.fcc((k3), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (256*sw2**2*lt.fcc((k1), (k2))*lt.fcc((k3), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (32*lt.fcc((k1), (k2))*lt.fcc((p1), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (128*sw2*lt.fcc((k1), (k2))*lt.fcc((p1), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (256*sw2**2*lt.fcc((k1), (k2))*lt.fcc((p1), (p1))*lt.fcc((p2), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + 2*(lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)) - lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (p2)) - lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((p1), (p2)) + lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((p1), (p2)) + lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((p2), (p2)) - lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((p2), (p2)))*(-4/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (16*sw2)/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2) + 2*(lt.fcc((k1), (p2))*lt.fcc((k2), (p1))*lt.fcc((k3), (k3)) - lt.fcc((k1), (p1))*lt.fcc((k2), (p2))*lt.fcc((k3), (k3)) - lt.fcc((k1), (p2))*lt.fcc((k2), (k3))*lt.fcc((k3), (p1)) + lt.fcc((k1), (k3))*lt.fcc((k2), (p2))*lt.fcc((k3), (p1)) + lt.fcc((k1), (p1))*lt.fcc((k2), (k3))*lt.fcc((k3), (p2)) - lt.fcc((k1), (k3))*lt.fcc((k2), (p1))*lt.fcc((k3), (p2)))*(4/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (16*sw2)/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + 4/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (16*sw2)/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + 4/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (16*sw2)/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2))))) + (2*lt.fcc((k1), (p2))*lt.fcc((k2), (p1)) - 2*lt.fcc((k1), (p1))*lt.fcc((k2), (p2)))*((-12*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (48*sw2*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (4*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (16*sw2*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (8*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (32*sw2*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (12*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (48*sw2*lt.fcc((k3), (k3)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (4*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (16*sw2*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (8*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (32*sw2*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (12*lt.fcc((k3), (k3)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (48*sw2*lt.fcc((k3), (k3)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (16*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (64*sw2*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (16*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (64*sw2*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (32*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (128*sw2*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2))))) + (2*lt.fcc((k1), (p1))*lt.fcc((k2), (k3)) - 2*lt.fcc((k1), (k3))*lt.fcc((k2), (p1)))*((-4*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (16*sw2*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (4*lt.fcc((p1), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (16*sw2*lt.fcc((p1), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (20*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (80*sw2*lt.fcc((k3), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (20*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (80*sw2*lt.fcc((p2), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (4*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (16*sw2*lt.fcc((k3), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (16*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (64*sw2*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2))))) + (2*lt.fcc((k1), (p2))*lt.fcc((k2), (k3)) - 2*lt.fcc((k1), (k3))*lt.fcc((k2), (p2)))*((20*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (80*sw2*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 - (20*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (80*sw2*lt.fcc((p1), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))**2 + (4*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (16*sw2*lt.fcc((k3), (p1)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 - (4*lt.fcc((p1), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (16*sw2*lt.fcc((p1), (p2)))/(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))**2 + (4*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (16*sw2*lt.fcc((k3), (p1)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) - (16*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))) + (64*sw2*lt.fcc((p1), (p2)))/((lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p1)) + lt.fcc((p1), (p1)))*(lt.fcc((k3), (k3)) - 2*lt.fcc((k3), (p2)) + lt.fcc((p2), (p2)))))))/(64*cw2**2*mz2**2*sw2**2)

    return a

def eeeegamma(vectors,masses=[0,0,0],s_sqrt = 1):
    # p = np.vstack((np.array([[0,0,0,0],[0.5,0,0.499999,0],[0.5,0,-0.499999,0]]),vectors))
    p = np.vstack((lt.four_vector_gen(E=s_sqrt,mass=0.005),vectors))
    p1 = lt.Lorentz4vector(components=p[1])
    p2 = lt.Lorentz4vector(components=p[2])
    p3 = lt.Lorentz4vector(components=p[3])
    p4 = lt.Lorentz4vector(components=p[4])
    p5 = lt.Lorentz4vector(components=p[5])

    coefficient = (4*np.pi*alpha)**3
    # s = lt.fcc((p[1]+p[2]),(p[1]+p[2]))
    s = (p1+p2)*(p1+p2)
    # s_ = lt.fcc((p[3]+p[4]),(p[3]+p[4]))
    s_ = (p3+p4)*(p3+p4)

    t13 = (p1-p3)*(p1-p3)
    t24 = (p2-p4)*(p2-p4)
    t14 = (p1-p4)*(p1-p4)
    t23 = (p2-p3)*(p2-p3)
    # print(t13)
    # print('t13')
    # print(t24)
    # print('t24')
    # print(s)
    # print('s')
    # print(s_)
    # print('s_')

    T = (1/(s*s_*t13*t24))*(s*s_*(s**2+s_**2)+t13*t24*(t13**2+t24**2)+t14*t23*(t14**2+t23**2))
    S = s/((p1*p5)*(p2*p5))+s_/((p3*p5)*(p4*p5))-t13/((p1*p5)*(p3*p5))-t24/((p2*p5)*(p4*p5))+t14/((p1*p5)*(p4*p5))+t23/((p2*p5)*(p3*p5))
    # print(T)
    # print('T')
    # print(S)
    # print('S')
    # print(lt.fcc(p[1],p[2]))
    # print(s)
    # print(s_)
    # print(coefficient*T*S)
    return [coefficient*T*S,p1,p2,p3,p4,p5,T,S]

# v = np.random.rand(3,4)
# eeeegamma(v)


def eemumu(vectors,masses=[0,0]):
    s_sqrt = 500
    E = s_sqrt/2
    s = s_sqrt**2
    p = np.vstack((np.array([[0,0,0,0],[250,0,0,250],[250,0,0,-250]]),vectors))
    coefficient = 2*np.pi*8*(e**4)/s_sqrt**4
    costheta = lt.cos_theta(p[1],p[3])
    numerator = 2*(lt.fcc(p[1],p[4])*lt.fcc(p[2],p[3])+lt.fcc(p[1],p[3])*lt.fcc(p[2],p[4])+
                 masses[1]**2*lt.fcc(p[1],p[2])+masses[1]**2*lt.fcc(p[3],p[4])+2*(masses[0]**2)*(masses[1]**2))
    # print(costheta)
    # print((2)/(137**2*4*250000))
    # t = masses[0]**2 + masses[1]**2-2*lt.fcc(p[2],p[4])
    # u = masses[0] ** 2 + masses[1] ** 2 - 2 * lt.fcc(p[2], p[3])
    # print(((s_sqrt/2)**4+(lt.v3_ip(p[1],p[3]))**2)/(137**2*(s_sqrt/2)**6*16))
    # p1 = lt.Lorentz4vector(components=p[1],mass=0)
    # p3 = lt.Lorentz4vector(components=p[3],mass=0)
    # numerator_1 = t ** 2 + u ** 2
    # print(numerator_1*coefficient*weight_all_zero/4)
    # print(numerator_1*(e**4)/(32*np.pi**2*s**3))
    # expression_1 = ((np.pi)**2*(8*E**4+8*lt.v3_ip(p[1],p[3])**2))/(137**2*2*np.pi**2*s**3)
    # # print(expression_1)
    # expression_2 = (E**4+lt.v3_ip(p[1],p[3])**2)/(137**2*16*(s_sqrt/2)**6)
    # print(expression_2)
    # expression_3 = (E**4+((p1.p3norm)*(p3.p3norm)*lt.cos_theta(p1,p3))**2)/(137**2*16*(s_sqrt/2)**6)
    # print(expression_3)
    # expression_4 = (E**4+(E*(p3.p3norm)*lt.cos_theta(p1,p3))**2)/(137**2*16*(s_sqrt/2)**6)
    # print(expression_4)
    # expression_5 = (1+lt.cos_theta(p1,p3)**2)/(137**2*4*s)
    # print(expression_5)
    # print(numerator/numerator_1)
    # print(p[1])
    # print(p[4])
    # print(coefficient*numerator*weight_all_zero)
    # print(((costheta**2+1)/(137**2*4*250000))/(coefficient*numerator*weight_all_zero))
    # print(p3.get_on_shell_ness())
    # print(p[3][3]**2+p[3][1]**2+p[3][2]**2)
    # print(lt.v3_ip(p[3],p[3]))
    return coefficient*numerator


# p1 = [0.5,0,0.499999,0]
# p2 = [0.5,0,-0.499999,0]
# p3 = [0.178093,-0.127916,0.050068,-0.113348]
# p4 = [0.356394,-0.028605,-0.183214,-0.304353]
# p5 = [0.465511,0.156521,0.133146,0.417701]
# v = np.array(p3)
# v = np.vstack((p3,p4))
# v = np.vstack((v,p5))
# print(eeeegamma(vectors=v)[0])


# print('!')
# print('bb')
# print(mc.phase_space_integration(eemumu,number_of_dots=10000,s_sqrt=500,masses=[0,0],dimension=2))
print(mc.phase_space_integration(eenunugamma1,number_of_dots=100000,s_sqrt=500,masses=[0,0,0],dimension=3,mode=0))
# vv = np.array([[45.5287,-16.429,7.88978,-41.7217],[249.981,81.8115,-45.7054,231.751],[204.49,-65.3826,37.8156,-190.029]])
# print(vv.shape)
# print(eenunugamma1(vectors=vv,masses=[0,0,0]))

