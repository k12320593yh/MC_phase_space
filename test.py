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
print(e)
gz = 0.718
gw = 0.629


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
    p = np.vstack((np.array([[0,0,0,0],[250,0,0,250],[250,0,0,-250]]),vectors))
    print(p)
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
    print(a)
    print(b)
    return coefficient*a/b

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





print(mc.phase_space_integration(eemumu,number_of_dots=10000,masses=[0,0,0],dimension=2))