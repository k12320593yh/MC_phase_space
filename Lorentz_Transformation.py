import numpy as np
import LorentzGenerators as lg
import scipy.linalg as sp
import sympy as sym

#Is it possible to define Lorentz Objects as ndarray's subclass?

#A general Lorentz 4 vector. Just what you learned in peskin.
class LorentzObject(object):
    def __init__(self,components,indices,cases):
        self.order = len(indices)
        self.components = np.array(components)
        self.indices = indices
        self.cases = cases
        # if self.components.shape == (4,4):
        #     self.ismetrictensor = True
        #     for i in range(4):
        #         for j in range(4):
        #             if self.components[i][j] != lg.minkowski_metric[i][j]:
        #                 self.ismetrictensor = False
        self.shape = self.components.shape
    # If we assume all indices are properly cased, then an 1D list would safely encode all necessary info
    # to describe the indices, whose elements are the indices name or the Lorentz Tensor and indices are the position
    # of each Lorentz index.
    def __mul__(self, other):
        # print("bing!")
        if not isinstance(other,LorentzObject):
            return LorentzObject(components=self.components*other,indices=self.indices,cases=self.case)
        else:
            return contract(self, other)

# class LorentzScalar(LorentzObject):
#     def __int__(self,value):
#         self.value = value
#
#     def __mul__(self, other):
#         if isinstance(other,LorentzScalar):
#             return self.value*other.value
#         else:
#             return self.value*other
#
#     def __add__(self, other):
#         if isinstance(other,LorentzScalar):
#             return self.value+other.value
#         else:
#             return self.value+other


#It is recommended to initialise a metric tensor with this to avoid complication.
class Metric_Tensor(LorentzObject):
    def __init__(self,indices):
        LorentzObject.__init__(self,components=lg.minkowski_metric,indices=['mu','nu'],cases=[True,True])
        self.ismetrictensor = True

    def __repr__(self):
        return lg.minkowski_metric

#This involves initialise and manipulate numpy ndarray with unknown shape
#May simply not be the correct path. What about turn to symbolic, algebraic solution?

class CaseError(Exception):
    def __int__(self):
        pass


def contract(tensor1, tensor2):
    new_indices = []
    tensor1_explicit = tensor1.components
    tensor2_explicit = tensor2.components
    tensor2_tmp = tensor2.components
    tensor1_contract_axes = []
    tensor2_contract_axes = []
    # print(tensor1.indices)
    # print(tensor2.indices)
    # Acquiring necessary axes information for contraction
    # new_indices keeps the indices of the result tensor
    # while the two contract_axes list keeps the ndarray axis index of the axes to be contracted over.
    for i in range(len(tensor1.indices)):
        if tensor1.indices[i] not in tensor2.indices:
            new_indices.append(i)
        else:
            j = tensor2.indices.index(tensor1.indices[i])
            if tensor1.cases[i] != tensor2.cases[j]:
                tensor1_contract_axes.append(i)
                tensor2_contract_axes.append(j)
            else:
                raise TypeError("Pair of to be contracted indices share common case.")
    for j in tensor2.indices:
        if j not in tensor1.indices:
            new_indices.append(j)
    if len(tensor2_contract_axes) == 0:
        return LorentzObject(components=np.tensordot(tensor1_explicit,tensor2_explicit,axes=0),indices=tensor1.indices+tensor2.indices,
                                 cases=tensor1.case+tensor2.cases)
    # If tensor 1 and tensor2 share no common indices, return the new tensor with
    # Lower the second tensor index before doing the contraction.
    # This might be useful for vague uses, but so far I haven't found a proper way to deal with metric
    # tensors in vague mode.
    """
    for i in tensor2_contract_axes:
        tensor2_tmp = np.tensordot(tensor2_tmp,lg.minkowski_metric,axes=([i],[0]))
    """
    # Now that all indices to be contracted over in tensor2 are placed in proper case,
    # do a normal tensor dot to get the result
    result_explict = np.tensordot(tensor1_explicit,tensor2_tmp,axes=(tensor1_contract_axes,tensor2_contract_axes))
    if len(new_indices) != 0:
        return LorentzObject(components=result_explict,indices=new_indices,cases=[True]*len(new_indices))
    else:
        return result_explict

# This function has an obvious flaw: it cannot handle the metric tensors.
# For now, all Lorentz objects are assumed to be wrongly upper cased and hence require
# a metric tensor to lower the index to have correct index case
#
# test_tensor_1 = LorentzObject(components=np.random.rand(4,4,4),indices=['a','b','c'])
# test_tensor_2 = LorentzObject(components=np.random.rand(4,4,4),indices=['e','f','g'])


#They're order-3 tensors while being order-1 Lorentz objects. Worth some caution.
class GammaMatricesObject(LorentzObject):
    def __init__(self,indices=['mu'],cases = [True],components=lg.gamma_matrices):
        if cases[0]:
            LorentzObject.__init__(self,components=components,indices=indices,cases=cases)
        else:
            LorentzObject.__init__(self,components=np.tensordot(lg.minkowski_metric,components,
                                                                axes=([1],[0])),indices=indices,cases=cases)

        # The last two dimensions are always used to store gamma indices.

    def __mul__(self, other):
        # print("ding!")
        return contractgamma(self,other)

class Gamma5(object):
    def __init__(self):
        pass
# g = GammaMatricesObject()
# print(g.fullcomponent)

def contractgamma(tensor1, tensor2):
    if isinstance(tensor2,Gamma5):
        tensor1_explicit = tensor1.components.reshape(4**len(tensor1.indices),4,4)
        for i in range(tensor1_explicit.shape[0]):
            tensor1_explicit[i] = np.matmul(tensor1_explicit[i],lg.gamma_5)
        result = np.reshape(tensor1_explicit,newshape=(4,)*len(tensor1.indices)+(4,4))
        return GammaMatricesObject(components=result,indices=tensor1.indices,cases=tensor1.cases)
    new_indices = []
    tensor1_explicit = tensor1.components
    tensor2_explicit = tensor2.components
    tensor2_tmp = tensor2.components
    tensor1_contract_axes = []
    tensor2_contract_axes = []
    # print(tensor1.indices)
    # print(tensor2.indices)
    # Acquiring necessary axes information for contraction
    # new_indices keeps the indices of the result tensor
    # while the two contract_axes list keeps the ndarray axis index of the axes to be contracted over.
    # This only goes through Lorentz indices so there is no need to rewrite.
    for i in range(len(tensor1.indices)):
        if tensor1.indices[i] not in tensor2.indices:
            new_indices.append(tensor1.indices[i])
        else:
            j = tensor2.indices.index(tensor1.indices[i])
            if tensor1.cases[i] != tensor2.cases[j]:
                tensor1_contract_axes.append(i)
                tensor2_contract_axes.append(j)
            else:
                raise TypeError("Pair of to be contracted indices share common case.")
    # tensor2_remaining_indices = len(tensor2.indices)-len(tensor2_contract_axes)
    for j in tensor2.indices:
        if j not in tensor1.indices:
            new_indices.append(j)
    if len(tensor2_contract_axes) == 0:
        # print((np.tensordot(tensor1_explicit,tensor2_explicit,axes=([-1],[-2]))).shape)
        # c = np.tensordot(tensor1_explicit, tensor2_explicit, axes=([-1], [-2])).swapaxes(-2, -3)
        c1 =np.tensordot(tensor1_explicit,tensor2_explicit,axes=([-1],[-2])).swapaxes(-2,-3)
        # print('dong!')
        # print(c[0,0])
        # print(c1[0,0])
        return GammaMatricesObject(components=c1,
                                   indices=tensor1.indices+tensor2.indices,
                                 cases=tensor1.cases+tensor2.cases)
    # If tensor 1 and tensor2 share no common indices, return the new tensor with
    # Lower the second tensor index before doing the contraction.
    # This might be useful for vague uses, but so far I haven't found a proper way to deal with metric
    # tensors in vague mode.
    """
    for i in tensor2_contract_axes:
        tensor2_tmp = np.tensordot(tensor2_tmp,lg.minkowski_metric,axes=([i],[0]))
    """
    # Now that all indices to be contracted over in tensor2 are placed in proper case,
    # do a normal tensor dot to get the result
    result_explict = np.tensordot(tensor1_explicit,tensor2_tmp,
                                  axes=(tensor1_contract_axes+[-1],tensor2_contract_axes+[-2]))
    if len(new_indices) != 0:
        return GammaMatricesObject(components=result_explict.swapaxes(-2,-3),indices=new_indices,cases=[True]*len(new_indices))
    else:
        return result_explict



class Lorentz4vector(LorentzObject):
    def __init__(self,components,mass = None,name = 'nothing',index = 'a',cases = True):
        LorentzObject.__init__(self,indices=index,components=components,cases=[cases])
        self.__name = name
#4-vector components
        self.components = np.array(components)
        self.index_name = index
#Particle mass. if you would like to operator at high energy limit, just set this to zero.
#If mass is not given, code will calculate it for you.
#Since we're talking about initial or final state particles mainly, off shell particles are not allowed.
        if mass == None:
            self.__mass = np.sqrt(four_vector_contraction_component(components,components))
#Designate the mass allows off shell particle but not recommended.
        else:
            self.__mass = mass
#Check if this 4-vector is on shell.
        self.__on_shell = check_on_shell_ness(self.components,self.__mass)
        self.transformation = lg.lorentz_generators[0]
#The four vectors are initialised as covariant 4 vectors.
        self.__original_components = components
        self.__mag = np.sqrt(self.components[0]**2-fcc(self.components,self.components))
        self.p3 = self.components[1:]
        self.e = self.components[0]
        self.shape = self.components.shape
        a = 0
        for i in range(3):
            a+=self.p3[i]**2
        self.p3norm = np.sqrt(a)

    # This is disabled in avoidance of confusion.
    # de__repr__(self):
    #     return self.componentsf

    def __str__(self):
        return str(self.components)

    def get_4_vector(self):
        return np.array(self.components)

    def __add__(self, other):
        if isinstance(other,Lorentz4vector):
            return Lorentz4vector(components=self.components+other.get_4_vector())
        else:
            return np.array(self.components+other)

    def __sub__(self, other):
        return self+other*(-1)

    def __mul__(self, other):
        if isinstance(other,Lorentz4vector):
            if self.indices == other.indices:
                return fcc(self.components,other.components)
            else:
                tmp = np.zeros(shape=(4,4))
                for i in range(4):
                    for j in range(4):
                        tmp[i][j] = self.components[i]*other.components[j]
                return LorentzObject(components=tmp,indices=list(self.indices+other.indices))
        elif isinstance(other,np.ndarray):
            return fcc(self.components,other)
        elif isinstance(other,GammaMatricesObject):
            return Slashed_Momentum(self.components)
        else:
            return Lorentz4vector(components=other*self.components)

    def sh(self,message=''):
        print(message+': '+str(self.components))

    def get_mass(self):
        return self.__mass

    def get_on_shell_ness(self):
        return self.__on_shell

    #Two modes: mode 0 will alter the original l4v and mode 1 will retun a new l4v object.
    def lorentz_transformation(self,ltf,mode=0):
        # tmp = np.matmul(lg.minkowski_metric,ltf.get_matrix_form())
        if mode == 0:
            self.components = np.matmul(ltf.get_matrix_form(),self.components)
        elif mode == 1:
            return Lorentz4vector(components=np.matmul(ltf.get_matrix_form(),self.components),mass=self.__mass)
        # self.__index_case = not self.__index_case

    def lorentz_transformation_off_matrix(self,matrix,mode=0):
        #Works exactly like the previous one but takes a np.ndarry as argument.
        if mode == 0:
            self.components = np.matmul(matrix,self.components)
        elif mode == 1:
            return Lorentz4vector(components=np.matmul(matrix,self.components),mass=self.__mass)

    #SO(3) rotate the momentum into z axis.
    def lorentz_rot_toz(self,record=False):
        rot_mat = find_toz_rotation(self.components)
        # tmp = np.matmul(lg.minkowski_metric,rot_mat)
        self.components = np.matmul(rot_mat,self.components)
        # self.transformation = np.matmul(rot_mat,self.transformation)
        # self.__index_case = not self.__index_case
        self.rotmat = rot_mat


    def zrapidity(self):
        if abs(self.components[1]) > 1e-8 or abs(self.components[2]) > 1e-8:
            raise TypeError("Three momentum not in z direction.")
        self.__zrapidity = find_boost_to_rest(self)
        if self.components[-1]<0:
            self.__zrapidity = -self.__zrapidity

    #Swap the case of index. p^\mu -> p_\mu, or vice versa.
    def swap_case(self):
        self.components = np.matmul(lg.minkowski_metric,self.components)
        self.index_case = not self.index_case

    def get_zrapidity(self):
        return self.__zrapidity

    #Must be done after rotate into z direction and
    def z_boost_to_rest_frame(self):
        boost_mat = boost_matrix(rapidities=np.array([0,0,-self.__zrapidity]))
        # print(boost_mat)
        # print(simple_z_boost(-self.__zrapidity))
        tmp = np.matmul(lg.minkowski_metric,boost_mat)
        self.components = np.matmul(tmp,self.components)
        # self.transformation = np.matmul(boost_mat,self.transformation)
        # self.__index_case = not self.__index_case
        self.boostmat = boost_mat


    def get_mag(self):
        return self.__mag

    def slashed(self):
        return Slashed_Momentum(self.components)





class Slashed_Momentum(object):
    def __init__(self,momentum):
        if isinstance(momentum,Lorentz4vector):
            p = momentum.get_4_vector()
        else:
            p = momentum
        self.matrix = p[0]*lg.gamma_0 - np.sum(p[i]*lg.gamma_matrices[i] for i in range(1,4))



def v3_ip(l4v1,l4v2):
    a = 0
    if isinstance(l4v1,Lorentz4vector):
        p1 = l4v1.get_4_vector()[1:]
    else:
        p1 = l4v1[1:]
    if isinstance(l4v2,Lorentz4vector):
        p2 = l4v2.get_4_vector()[1:]
    else:
        p2 = l4v2[1:]

    for i in range(3):
        a += p1[i]*p2[i]
    return a

def cos_theta(l4v1,l4v2):
    if isinstance(l4v1,Lorentz4vector):
        t_componets_1 = l4v1.get_4_vector()[1:]
    else:
        t_componets_1 = l4v1[1:]
    if isinstance(l4v2,Lorentz4vector):
        t_componets_2 = l4v2.get_4_vector()[1:]
    else:
        t_componets_2 = l4v2[1:]
    a = 0
    b = 0
    c = 0
    for i in range(3):
        a += t_componets_1[i]*t_componets_2[i]
        b += t_componets_1[i]*t_componets_1[i]
        c += t_componets_2[i]*t_componets_2[i]
    return a/np.sqrt(b*c)


def four_vector_add(l4v1,l4v2):
    summ = l4v1.get_4_vector()+l4v2.get_4_vector()
    #p^2 = m^2 for every on shell particle
    pmass = np.sqrt(summ[0]**2 - np.sum(summ[i]**2 for i in range(1,4)))
    return Lorentz4vector(components=summ,name='new vector',mass=pmass)


#4vector contracion. Returning a lorentz invariant scalar.
#Given that most people will not initialise Lorentz indices case with caution
#The requirement for cases of the to-be-contracted over index to be different is removed.
def four_vector_contraction(l4v1,l4v2):
    if (l4v1.index_name == l4v2.index_name):
        v1 = l4v1.get_4_vector()
        v2 = l4v2.get_4_vector()
        return four_vector_contraction_component(v1,v2)


def four_vector_contraction_component(l4v1,l4v2):
    return l4v1[0]*l4v2[0]-np.sum(l4v1[i]*l4v2[i] for i in range(1,4))


def three_vector_gen_off_t_p(p,theta,phi):
    return [np.sin(theta)*np.cos(phi)*p,np.sin(theta)*np.sin(phi)*p,np.cos(theta)*p]




def fcc(l4v1,l4v2):
    return l4v1[0] * l4v2[0] - np.sum(l4v1[i] * l4v2[i] for i in range(1, 4))


def check_on_shell_ness(components,mass):
    a = components[0] ** 2 - components[1] ** 2 - components[2] ** 2 - components[3] ** 2 - mass ** 2
    if a <= 1e-8:
        return True
    else:
        return False


#A slight differece from zero can be viewed as error generated in numerical calculations.
#IS 1e-8 too loose? Not so sure.


#Generate Lorentz group elements using generators $J_i$ and $K_i$
#$\theta_i$ and $\beta_i$ are required as input
def general_matrix_form(parameter):
    summ = np.zeros((4,4))
    for i in range(1,7):
        summ -= parameter[i-1]*lg.lorentz_generators[i]
        matrix = sp.expm(summ)
    return np.array(matrix)


#Calculating |p|
def spatial_norm(vector):
    norm = np.sqrt(np.sum(i ** 2 for i in vector[1:]))
    return norm


def rot_mat_xy(phi):
    rot_mat = np.diag((1.,1.,1.,1.))
    rot_mat[1,1] = np.cos(phi)
    rot_mat[2,2] = np.cos(phi)
    rot_mat[1,2] = np.sin(phi)
    rot_mat[2,1] = -np.sin(phi)
    return rot_mat

def rot_mat_xz(phi):
    rot_mat = np.diag((1.,1.,1.,1.))
    rot_mat[1,1] = np.cos(phi)
    rot_mat[3,3] = np.cos(phi)
    rot_mat[1,3] = -np.sin(phi)
    rot_mat[3,1] = np.sin(phi)
    return rot_mat


# Find the Lorentz Transformation in matrix form that will transform the given vector into z direction
# with the first component intact.
def find_toz_rotation(vector):
    if vector[1] == vector[2] == 0:
        return lg.lorentz_generators[0]
    #Acquiring the rotation angles to z. Just like SO(3)
    norm = spatial_norm(vector)
    phi = np.arccos(vector[1]/np.sqrt(vector[1]**2+vector[2]**2))
    theta = np.arccos(vector[3]/norm)
    #Calculationg the rotation martix with naive matrices multiplication
    # mat_rotphi = sp.expm(-phi*lg.lorentz_generators[3])
    # mat_rottheta = sp.expm(-theta*lg.lorentz_generators[2])
    matrix = np.matmul(rot_mat_xz(theta),rot_mat_xy(phi))
    # print(matrix)
    return matrix

def boost_matrix(rapidities):
    summ = np.zeros((4,4))
    for i in range(3):
        summ -= lg.lorentz_generators[i+4]*rapidities[i]
    matrix = sp.expm(-summ)
    return matrix

def simple_z_boost(rapidity):
    mat = np.diag((1.,1.,1.,1.))
    mat[0,0] = mat[3,3] = np.cosh(rapidity)
    mat[0,3] = mat[3,0] = np.sinh(rapidity)
    return mat


def find_boost_to_rest(l4v):
    if not l4v.get_on_shell_ness():
        print('Warning: OFF SHELL PARTICLE IN FINAL STATE.')
    if not l4v.get_mass():
        print('IS MASSLESS PARTICLE')
        return np.Infinity
    vector = l4v.get_4_vector()
    try:
        rapidity = 0.5*np.log((vector[0]+spatial_norm(vector))/(vector[0]-spatial_norm(vector)))
    except ZeroDivisionError:
        print('Warning: OFF SHELL PARTICLE IN FINAL STATE.')
        return np.NaN
    # print(rapidity)
    return rapidity

def four_vector_gen(E,mass):
    axial_compoent = np.sqrt(E**2-mass**2)
    return np.array([[0,0,0,0],[E,0,0,axial_compoent],[E,0,0,-axial_compoent]])

#General Lorentz Transformations Generated by Two Sets of Generators
class LorentzTransFormation(object):
    def __init__(self,name,parameter=None,matrix=None):
        self.__name = name
        if type(parameter) != None:
            self.__parameter = parameter
            self.__matrix_form = general_matrix_form(self.__parameter)
        elif type(matrix) != None & type(parameter) == None:
            self.__matrix_form = matrix
        else:
            raise TypeError("No valid datatype to initialise a LorentzTransFormation object")

    def show_matrix_form(self):
        print(self.__matrix_form)

    def get_matrix_form(self):
        return self.__matrix_form

    #backdoor measure to overwrite the matrix of a LorentzTransFormation object.
    #If you're going to use this method, initial the stance with a straight 0 parameter ndarray.
    def set_matrix_form(self,matrix):
        self.__matrix_form = matrix


def gamma_matrices_trace_generator(gamma):
    if not isinstance(gamma,GammaMatricesObject):
        raise TypeError("Not a gamma object")
    if len(GammaMatricesObject.indices)%2 != 0:
        return 0
    else:
        if len(GammaMatricesObject.indices) == 2:
            return Metric_Tensor(indices=GammaMatricesObject)
        else:
            gamma_matrices_trace_generator()

# ga_mu = GammaMatrices(index='mu')
# p0 = Lorentz4vector(components=[4,3,2,1])
# ps = p0.slashed()
# print((p0*ga_mu).matrix)
# print(ps.matrix)
# print(np.trace(ps.matrix))
# print(p0*p0)
# a = Lorentz4vector(components=[4,3,2,1],index='a')
# b = Lorentz4vector(components=[10,4,4,4],index='b')
# print((a*b).components)

# p1 = Lorentz4vector(components=[1,0,0,1])
# p2 = Lorentz4vector(components=[2,0,0,1])

# F = 0
# for i in range(10000):
#     dots = np.random.rand(2,)
#     dots[0] *= np.pi
#     dots[1] *= 2*np.pi
#     ap = three_vector_gen_off_t_p(1,dots[0],dots[1])
#     a = Lorentz4vector(components=[1]+ap,mass=0)
#     if cos_theta(a,[1,0,0,1]) - np.cos(dots[0]) > 1e-6:
#         F += 1
# print(F)