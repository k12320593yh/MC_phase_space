import sympy as sym
import scipy.linalg as sp
import threeD_rotation as td
import numpy as np
import Lorentz_Transformation as lt
import LorentzGenerators as lg
import MonteCarlo as mc
import time as time

test_tensor_1 = lg.minkowski_metric
test_tensor_2 = np.random.rand(4,4,4)
test_tensor_3 = np.random.rand(4,4,4)
p1 = lt.LorentzObject(components=[5,4,3,2],indices=['mu'])
p2 = lt.LorentzObject(components=[6,5,3,4],indices=['mu'])
g1 = lt.LorentzObject(components=lg.minkowski_metric,indices=['mu','nu'])
g2 = lt.LorentzObject(components=lg.minkowski_metric,indices=['mu','nu'])
a = np.reshape(np.linspace(0,63,64),(4,4,4))
print(a)
# print(test_tensor_2[:,:,0])
# print(test_tensor_2.take(0,axis=2))
#np.take
c = np.tensordot(test_tensor_2,test_tensor_3,axes=([0],[0]))
print(c.shape)
print(lt.contract(p1,p2).components)

# print("Running time : %f CPU seconds" % (end_CPU - start_CPU))
