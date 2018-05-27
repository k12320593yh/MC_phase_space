import numpy as np
import Lorentz_Transformation as lt
import LorentzGenerators as lg

#This module is dedicated to do symbolic operations, manipulations, and calculations of Lorentz objects
#Without involving explicit form.

#m for mu, n for nu, a for alpha, b for beta
#For now this only takes care of Lorentz four momenta, gamma matrices and some ordinary constant tensors
#like levi-civita.


class SymbolicLorentzObject(object):
    def __init__(self,indices=['m'],cases=[True]):
        self.index = indices
        self.case = cases

    def __mul__(self, other):
        pass

class SFourVector(SymbolicLorentzObject):
    def __init__(self,index='m',case=True):
        SymbolicLorentzObject.__init__(self,indices=[index],case=[case])

    def __mul__(self, other):
        if self.index == other.index:
            return SymbolicLorentzObject(indices=self.index+other.index,)

    def explict(self,expilcit_l4v):
        if isinstance(expilcit_l4v,lt.Lorentz4vector):


class ScalarProduct(object):
    def __init__(self,v1,v2):
        self.v1 = v1
        self.v2 = v2

    def value(self):
        return self.v1.e*self.v2.e

p = SymbolicLorentzObject()
q = SymbolicLorentzObject()