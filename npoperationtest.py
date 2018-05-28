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
k1 = np.reshape(np.linspace(0,15,16),newshape=(4,4))
k2 = np.array([5,6,7,8])
# print(np.tensordot(k1,k2,axes=0).shape)
p1 = lt.LorentzObject(components=[5,4,3,2],indices=['mu'],cases=[True])
p2 = lt.LorentzObject(components=[6,5,3,4],indices=['mu'],cases=[False])
p3 = lt.LorentzObject(components=[7,-1,-1,1],indices=['mu'],cases=[False])
g1 = lt.LorentzObject(components=lg.minkowski_metric,indices=['mu','nu'],cases=[False,False])
g2 = lt.LorentzObject(components=lg.minkowski_metric,indices=['a','b'],cases=[True,True])
a = np.reshape(np.linspace(0,63,64),(4,4,4))
ga1 = lt.GammaMatricesObject(indices=['mu'])
ga2 = lt.GammaMatricesObject(indices=['nu'])
ga3 = lt.GammaMatricesObject(indices=['alpha'])
ga4 = lt.GammaMatricesObject(indices=['beta'])
ga5 = lt.Gamma5()
ga6 = lt.GammaMatricesObject(indices=['mu'],cases=[False])
print(ga1.components[...,0,3])
# print(ga1.components[...,1,1])
# print(ga2.components[1])
# print(a)
# print(test_tensor_2[:,:,0])
# print(test_tensor_2.take(0,axis=2))
#np.take
# c = np.tensordot(test_tensor_2,test_tensor_3,axes=([0],[0]))
# print(c.shape)
# d = ga1*ga2
# print(d.components[0,0])
# e = d*ga3
# print(e.components.shape)
# print(e.components[0,0,0])
# f = e*ga4
# print(f.components[1,1,0,0])
# print(f.indices)
# g = f*ga5
# print(g.components.shape)
# print(g.indices)
# print(g.components[0,0,0])

# print(d.components[0,1])
#So far this works fine for two gamma objects.
# print(np.dot(lg.gamma_0,lg.gamma_1))
# print("Running time : %f CPU seconds" % (end_CPU - start_CPU))
# print((ga1*ga2).indices)
# print(ga1.components.shape)
# print(ga2.components.shape)
# print((np.tensordot(ga1.components,ga2.components,axes=([-1],[-2]))).shape)
# ga1ga2 = ga1*ga2
# print(ga1ga2.components.shape)
# print(ga1ga2.components.shape)
# print(ga3.shape)
# print(isinstance(ga1ga2,lt.GammaMatricesObject))
# c_g12 = ga1ga2.components
# print(c_g12.shape)
# c_g3 = ga3.components
# print(c_g3.shape)
# print((np.tensordot(c_g12,c_g3,axes=([-1],[-2])).shape))
# print((ga1ga2*ga3).components[0,0,0])