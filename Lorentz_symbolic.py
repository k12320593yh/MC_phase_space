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
        self.coefficients = objects[0].coefficients
        self.names = [objects[0].names[:]]
        # print(objects)
        if len(objects) > 1:
            for i in range(1,len(objects)):
                for j in range(len(objects[i].indices)):
                    if isinstance(objects[i],NonLorentzObject):
                        for k in range(len(objects[i].coefficients)):
                            if objects[i].coefficients[k] != 1:
                                self.coefficients.append(objects[i].coefficients[k])
                            continue
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
                            # print(objects[i].names)
                            # print(objects[i].cases)
                            if isinstance(objects[i].cases,list):
                                self.cases.append(objects[i].cases[j])
                            else:
                                self.cases.append(objects[i].cases)
                            self.names.append(objects[i].names)
                            for k in range(len(objects[i].coefficients)):
                                if objects[i].coefficients[k] != 1:
                                    self.coefficients.append(objects[i].coefficients[k])

    # Returns a human readable string.
    def __str__(self):
        string = ''
        for i in range(len(self.objects)):
            string+=str(self.objects[i].names)
            # string+=' '
            for j in range(len(self.objects[i].indices)):
                if self.objects[i].cases[j]:
                    string+='^'
                    string+=str(self.objects[i].indices[j])
                else:
                    string+='_'
                    string+=str(self.objects[i].indices[j])
            string+=" "
        return str(string)


    # If is abelian object like gmunu,pmu, then search for object in the sequence to contract
    def __mul__(self, other):
        if isinstance(other,NonLorentzObject):
            self.coefficients.append(other)
            return self


        if other.isabelian:
            print(other.cases)
            # Is this an internal bug of python?
            # I did nothing that might change the value of other.cases,
            # however it was changed.
            tmp_objects = self.objects[:]
            for i in range(len(self.objects)):
                for j in range(len(self.objects[i].indices)):
                    if self.objects[i].indices[j] in other.indices:
                        print(type(tmp_objects[i]))
                        print(other.cases)
                        tmp_objects[i] = tmp_objects[i]*other
                        return LorentzSequence(tmp_objects)
    # If not, only try the last element.
        else:
            tmp_objects = self.objects[:]
            tmp_objects[-1] = tmp_objects[-1]*other
            return LorentzSequence(objects=tmp_objects+other.objects)

    # if isinstance(other,LorentzSequence):
        return LorentzSequence(objects=self.objects+other.objects)



class SymbolicLorentzObject(object):
    def __init__(self, indices, cases, names,coefficients=1):
        self.indices = indices
        self.cases = cases
        self.names = names
        self.objects = [self]
        self.coefficients = [coefficients]
        self.isabelian = True

class NonLorentzObject(object):
    def __init__(self, names, value):
        self.names = names
        self.value = value
        self.coefficient = value


class Scalar(NonLorentzObject):
    def __init__(self, value=1, names='const'):
        NonLorentzObject.__init__(self, value=value, names=names)


    def __mul__(self, other):
        return Scalar(value=(self.value*other.value))

    def __add__(self, other):
        return Scalar(value=(self.value+other.value))

    def __str__(self):
        return str(self.value)



class ScalarProduct(Scalar):
    def __init__(self,S4V1,S4V2):
        Scalar.__init__(self,value=1)
        self.members = [S4V1,S4V2]

class Vector(NonLorentzObject):
    def __init__(self,value,names):
        NonLorentzObject.__init__(self,value=value,names=names)

    def __mul__(self, other):
        return self.value*other


class SFourVector(SymbolicLorentzObject):
    def __init__(self,indices,cases,names='p',coefficients=1):
        # print(index)
        SymbolicLorentzObject.__init__(self,indices=indices, cases=cases, names=names, coefficients=coefficients)


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


        if isinstance(other,MetricTensor):
            # print('Ding!')
            return other*self


class GammaFourVector(SFourVector):
    #name of gamma four vector are forced to be gamma.
    def __init__(self,indices):
        SFourVector.__init__(self,indices=indices,cases=True,names='gamma')
        self.isabelian = False




def mtinitialiser(indices,cases):
    if indices[0] != indices[1]:
        # print('bing')
        return MetricTensor(indices=indices,cases=cases)
    else:
        # print('ding!')
        if cases[0] == cases[1]:
            return NonLorentzObject(value=-2,names='g')
        else:
            return NonLorentzObject(value=4,names='g')


class MetricTensor(SymbolicLorentzObject):
    def __init__(self,indices,cases):
        SymbolicLorentzObject.__init__(self,indices=indices,cases=cases,names='g',coefficients=1)

    def __mul__(self, other):
    #This part handles contranction with four vectors like pmu.
        if isinstance(other,SFourVector):
            for i in range(2):
                # print(other.indices[0])
                # print(self.indices[i])
                if self.indices[i] == other.indices[0]:
                    # print(other.cases[0])
                    if self.cases[i] != other.cases[0]:
                        new_indices = self.indices[:]
                        # print(new_indices)
                        new_indices.pop(i)
                        new_cases = self.cases[:]
                        new_cases.pop(i)

                        return SFourVector(indices=new_indices,cases=new_cases[0],names=other.names)
                    else:
                        raise lt.CaseError("Trying to lower a lower index or raise an upper index.")
            print("Ding!")
            return LorentzSequence(objects=[self,other])



        if isinstance(other,MetricTensor):
            for i in range(2):
                for j in range(2):
                    if self.indices[i] == other.indices[j]:
                        if self.cases[i] == other.cases[j]:
                            raise lt.CaseError("Trying to lower a lower index or raise an upper index.")
                        else:
                            tmp_indices_1 = self.indices[:]
                            tmp_indices_2 = other.indices[:]
                            tmp_indices_1.pop(i)
                            tmp_indices_2.pop(j)
                            new_indices = tmp_indices_1+tmp_indices_2
                            # print(new_indices)
                            tmp_cases_1 = self.cases[:]
                            tmp_cases_2 = other.cases[:]
                            tmp_cases_1.pop(i)
                            tmp_cases_2.pop(j)
                            new_cases = tmp_cases_1+tmp_cases_2
                            # print(new_cases)
                            return mtinitialiser(indices=new_indices,cases=new_cases)
            return LorentzSequence(objects=[self,other])
p1 = SFourVector(indices=['m'],names='p1',cases=[False])
g1 = mtinitialiser(indices=['a','n'],cases=[False,False])
g1_ = mtinitialiser(indices=['a','m'],cases=[True,True])
g2 = mtinitialiser(indices=['b','c'],cases=[True,True])
g3 = mtinitialiser(indices=['d','d'],cases=[False,True])
s2 = Scalar(2)
print((p1*g1_).indices)
# print(type(g3))
# print(g3.value)
# print(type(g1*g1_))
# print((g1*g1_).cases)
# print((g1*g1_).indices)
# print(type(g1_*g1))
x = (g2*g1)*g1_
print((x*p1).names)
# print(type((g2*g1)*g1_))
# print(x)
# print(x.cases)
# print(type(p1*g1_))
# print((p1*g1_).names)
# print((g1_*g1).indices)
# print((p1*g1).indices)

# print((p1*g1).cases)
# print((p1*g1).objects[0])
# s1 = Scalar(1)
# s2 = Scalar(2)
# print(s1*s2)
# print(p1.indices)
# p2 = SFourVector(indices=['n'],names='p2')
# p3 = SFourVector(indices=['a'],names='p3')
# p12 = p1*p2
# print(p12.names)
# print(p12.cases)
# print(p12.indices)
# p123 = p12*p3
# print(p123.cases)
# print(p123.indices)
# print(p123.names)
# p312 = p3*p12
# print(p312.indices)
# print(p312.names)
# p1*(p1*p2)
# print(p1)
# print([p1,p2])
# print((p1*p2).indices)
# print(p3.indices)
# # print((p1*p2))
# class MinkowskiMetric(SymbolicLorentzObject):
#     def __init__(self,indices=['mu','nu'],cases=[True,True]):
#         SymbolicLorentzObject.__init__(self,indices=indices,cases=cases)
#         pass
#
#
#
#
# class GammaObject(SymbolicLorentzObject):
#     def __init__(self,index='m',case=True):
#         SymbolicLorentzObject.__init__(self,indices=[index],case=[case])
#
# class ScalarProduct(object):
#     def __init__(self,v1,v2):
#         self.v1 = v1
#         self.v2 = v2
#
#     def value(self):
#         return self.v1.e*self.v2.e

# p = SymbolicLorentzObject()
# q = SymbolicLorentzObject()