import sympy as sp
import numpy as np
import LorentzGenerators as lg
import Lorentz_Transformation as lt
# np.random.seed(42)

s_sqrt,masses = 500,(0,0)
# E = 0
# for i in range(10000):
raw_dot = np.random.rand(5,)
mv = raw_dot[0]*(s_sqrt-np.sum(masses))+np.sum(masses)
E_gamma = (s_sqrt**2-mv**2)/(2*s_sqrt)
theta_gamma = np.pi*raw_dot[1]
phi_gamma = 2*np.pi*raw_dot[2]
print(lg.lorentz_generators[1:4])
# print(theta_gamma)
# print(phi_gamma)
# E += E_gamma
# print(E)
E_1 = raw_dot[3]*(s_sqrt-E_gamma)
# It's a nice approach to generate photon energy using chain decay idea.
# However in this case all we need is E_gamma.