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




# print("Running time : %f CPU seconds" % (end_CPU - start_CPU))
