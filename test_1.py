import LorentzGenerators as lg
import Lorentz_Transformation as lt
# import MonteCarlo as mc
import numpy as np
import scipy.linalg as sp
import inspect as insp

p1 = [1,0,0,0]
p2 = [2,0,0,0]
print(lt.fcc(p1,p2))