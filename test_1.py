import LorentzGenerators as lg
import Lorentz_Transformation as lt
# import MonteCarlo as mc
import numpy as np
import scipy.linalg as sp
import inspect as insp

x = np.linspace(0,80,81)
# print(x)
b = x.reshape(3,3,3,3)
a = np.random.rand(3,3,3)
print(np.tensordot(a,b,axes=([0],[0])))
# print(b)
print(b)
print("############")
# print(b[:,:,0])
print("########")
# print(b[0,...])
print(b.transpose(2,0,1,3)[0])
print("########")
# print(b.swapaxes(0,2)[0])