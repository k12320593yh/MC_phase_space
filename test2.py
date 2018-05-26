import sympy as sym
import scipy.linalg as sp
import threeD_rotation as td
import numpy as np
import Lorentz_Transformation as lt
import LorentzGenerators as lg
import MonteCarlo as mc
import time as time

test_decay0 = lt.Lorentz4vector(name='test',components=[5,0,0,3],mass=4)
test_decay1 = lt.Lorentz4vector(name='test',components=[5,0,0,-3],mass=4)
test_masses = [1,2,3]

start_CPU = time.clock()
for i in range(10):
    p1 = mc.two_three_phase_space_dot(momentum_0=test_decay0,momentum_1=test_decay1,masses=test_masses)
    print(np.sum(p1.get_vector()[2:,0]))
end_CPU = time.clock()
#Generate 10000events takes ~45s
#Not bad!



print("Running time : %f CPU seconds" % (end_CPU - start_CPU))
