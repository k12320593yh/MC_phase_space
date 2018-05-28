import numpy as np
import Lorentz_Transformation as lt
import LorentzGenerators as lg

#This module is dedicated to do symbolic operations, manipulations, and calculations of Lorentz objects
#Without involving explicit form.

#m for mu, n for nu, a for alpha, b for beta
#For now this only takes care of Lorentz four momenta, gamma matrices and some ordinary constant tensors
#like levi-civita.

#Lorentz Objects ADDED together.
#Designed to be a list comprising multiple LorentzSequences with same Lorentz indices.

class LorentzObjectComplex(object):
    def __init__(self,sequences):
        self.sequences = []
        self.required_indices = sequences[0].indices
        for i in range(len(sequences)):
            if sequences[i].indices != self.required_indices:
                raise IndexError("Trying to add object with different Lorentz indices together.")
            else:
                self.sequences.append(sequences[i])

#A Sequence of Lorentz object dotted together.
class LorentzSequence(object):
    def __init__(self,objects):
        self.indices = objects[0].indices[:]
        self.cases = objects[0].cases[:]
        self.objects = objects
        self.coefficient = [1,]
        self.names = []
        # print(objects)
        if len(objects) > 1:
            for i in range(1,len(objects)):
                for j in range(len(objects[i].indices)):
                    if isinstance(objects[i],Scalar):
                        self.coefficient.append(objects[i])
                    else:
                        if objects[i].indices[j] in self.indices:
                            if objects[i].cases[j] == self.cases[self.indices.index(objects[i].indices[j])]:
                                raise lt.CaseError("Trying to construct a Lorentz Object with duplicate indices.")
                            #Denies sequences like p^mu g^munu, allows p^mu g_munu
                            else:
                                print("Contraction attention")
                                pass
                            #will take care of this later.

                        else:
                            self.indices.append(objects[i].indices[j])
                            self.cases.append(objects[i].cases[j])
                            self.names.append(objects[i].names)

    def __str__(self):
        string = ''
        for i in range(len(self.objects)):
            string+=self.objects[i].names
            string+="^"
            for j in range(len(self.objects[i].indices)):
                string+=self.objects[i].indices[j]
            string+=" x "
        return str(string)

class SymbolicLorentzObject(object):
    def __init__(self,indices=['m'],cases=[True],name='p'):
        self.indices = indices
        self.cases = cases
        self.names = name

class Scalar(object):
    def __init__(self,value = 1):
        self.value = value


class ScalarProduct(object):
    def __init__(self,S4V1,S4V2):
        self.members = [S4V1,S4V2]

class SFourVector(SymbolicLorentzObject):
    def __init__(self,indices=['m'],case=True,names='p'):
        # print(index)
        SymbolicLorentzObject.__init__(self,indices=indices,cases=[case],name='p')

    def __str__(self):
        return self.names

    def __mul__(self, other):
        if isinstance(other,SFourVector):
            if self.indices[0] != other.indices[0]:
                return LorentzSequence(objects=[self,other])
            elif self.indices[0] == other.indices[0]:
                if self.cases[0] == other.cases[0]:
                    return ScalarProduct(self,other)
                else:
                    raise lt.CaseError("Trying to construct a Lorentz Object with duplicate indices.")

        if isinstance(other,LorentzSequence):
            return LorentzSequence(objects=[self]+other.objects)

p1 = SFourVector(indices=['m'])
print(p1.indices)
p2 = SFourVector(indices=['n'])
# p3 = SFourVector(indices='a')
# print(p1)
# print([p1,p2])
# print((p1*p2).indices)
# print(p3.indices)
print((p1*p2))
class MinkowskiMetric(SymbolicLorentzObject):
    def __init__(self,indices=['mu','nu'],cases=[True,True]):
        SymbolicLorentzObject.__init__(self,indices=indices,cases=cases)
        pass




class GammaObject(SymbolicLorentzObject):
    def __init__(self,index='m',case=True):
        SymbolicLorentzObject.__init__(self,indices=[index],case=[case])

class ScalarProduct(object):
    def __init__(self,v1,v2):
        self.v1 = v1
        self.v2 = v2

    def value(self):
        return self.v1.e*self.v2.e

p = SymbolicLorentzObject()
q = SymbolicLorentzObject()