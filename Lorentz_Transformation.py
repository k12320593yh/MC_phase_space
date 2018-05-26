import numpy as np
import LorentzGenerators as lg
import scipy.linalg as sp

#A general Lorentz 4 vector. Just what you learned in peskin.
class LorentzObject(object):
    def __init__(self,order):
        self.order = order


class Lorentz4vector(LorentzObject):
    def __init__(self,components,mass = None,name = 'nothing',index_name = 'mu'):
        LorentzObject.__init__(self,order=1)
        self.__name = name
#4-vector components
        self.__components = np.float64(components)
#Particle mass. if you would like to operator at high energy limit, just set this to zero.
#If mass is not given, code will calculate it for you.
#Since we're talking about initial or final state particles mainly, off shell particles are not allowed.
        if mass == None:
            self.__mass = np.sqrt(four_vector_contraction_component(components,components))
#Designate the mass allows off shell particle but not recommended.
        else:
            self.__mass = mass
#Check if this 4-vector is on shell.
        self.__on_shell = check_on_shell_ness(self.__components,self.__mass)
        self.transformation = lg.lorentz_generators[0]
#The four vectors are initialised as covariant 4 vectors.
        self.index_case = [True]
        self.__original_components = components
        self.index_name = [index_name]
        self.__mag = np.sqrt(self.__components[0]**2-fcc(self.__components,self.__components))
        self.p3 = self.__components[1:]
        self.e = self.__components[0]
        a = 0
        for i in range(3):
            a+=self.p3[i]**2
        self.p3norm = np.sqrt(a)
    def get_4_vector(self):
        return np.array(self.__components)

    def sh(self,message=''):
        print(message+': '+str(self.__components))

    def get_mass(self):
        return self.__mass

    def get_on_shell_ness(self):
        return self.__on_shell

    #Two modes: mode 0 will alter the original l4v and mode 1 will retun a new l4v object.
    def lorentz_transformation(self,ltf,mode=0):
        # tmp = np.matmul(lg.minkowski_metric,ltf.get_matrix_form())
        if mode == 0:
            self.__components = np.matmul(ltf.get_matrix_form(),self.__components)
        elif mode == 1:
            return Lorentz4vector(components=np.matmul(ltf.get_matrix_form(),self.__components),mass=self.__mass)
        # self.__index_case = not self.__index_case

    def lorentz_transformation_off_matrix(self,matrix,mode=0):
        #Works exactly like the previous one but takes a np.ndarry as argument.
        if mode == 0:
            self.__components = np.matmul(matrix,self.__components)
        elif mode == 1:
            return Lorentz4vector(components=np.matmul(matrix,self.__components),mass=self.__mass)

    #SO(3) rotate the momentum into z axis.
    def lorentz_rot_toz(self,record=False):
        rot_mat = find_toz_rotation(self.__components)
        # tmp = np.matmul(lg.minkowski_metric,rot_mat)
        self.__components = np.matmul(rot_mat,self.__components)
        # self.transformation = np.matmul(rot_mat,self.transformation)
        # self.__index_case = not self.__index_case
        self.rotmat = rot_mat


    def zrapidity(self):
        if abs(self.__components[1]) > 1e-8 or abs(self.__components[2]) > 1e-8:
            raise TypeError("Three momentum not in z direction.")
        self.__zrapidity = find_boost_to_rest(self)
        if self.__components[-1]<0:
            self.__zrapidity = -self.__zrapidity

    #Swap the case of index. p^\mu -> p_\mu, or vice versa.
    def swap_case(self):
        self.__components = np.matmul(lg.minkowski_metric,self.__components)
        self.index_case = not self.index_case

    def get_zrapidity(self):
        return self.__zrapidity

    #Must be done after rotate into z direction and
    def z_boost_to_rest_frame(self):
        boost_mat = boost_matrix(rapidities=np.array([0,0,-self.__zrapidity]))
        # print(boost_mat)
        # print(simple_z_boost(-self.__zrapidity))
        tmp = np.matmul(lg.minkowski_metric,boost_mat)
        self.__components = np.matmul(tmp,self.__components)
        # self.transformation = np.matmul(boost_mat,self.transformation)
        # self.__index_case = not self.__index_case
        self.boostmat = boost_mat


    def get_mag(self):
        return self.__mag


    # def get_indice(self):
    #     pass


class Lorentz_order_2_tensor(LorentzObject):
    def __init__(self,components = lg.minkowski_metric,names=['mu','nu'],cases=[True,True],name = 'order 2 Lorentz Tensor'):
        LorentzObject.__init__(self,order=2)
        self.components = components
        self.names = names
        self.cases = cases
        indices = []
        for i in range(2):
            indices.append((names[i], [i, cases[i]]))
        self.indices = dict(indices)

    def expand_along_axis(self,axis):
        if axis == 0:
            return self.components
        elif axis == 1:
            return self.components.transpose()


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


#Not Ready yet. Do not use.
class Lorentz_HD_tensor(object):
#Initialised to be g^{\mu\nu} by default.
#Two attributes are initialised with the object: An ndarray storing explict matrix form, another is a dict
#whose keys are the index names and values binary tuples, 1st component being the location while the second being
#the case.
    def __init__(self,name = 'HDtensor',components = lg.minkowski_metric,names = ['mu','nu'],cases = [True,True]):
        self.components = components
        if (len(names) != len(cases)) or len(names) != len(components.shape) or len(cases) != len(components.shape):
            raise TypeError("Tensor dimensionality does not match number of indices.")
        indices = []
        for i in range(len(names)):
            indices.append((names[i],[i,cases[i]]))
        self.indices = dict(indices)


    def get_indice(self):
        return self.indices

    def expand_to_component_vectors_along_certain_index(self,index):
        pass
        # axis_number = self.indices[index][0]
        # components_tmp = np.swapaxes(self.components,0,axis_number)
        # vectors = []
        # for i in range(components_tmp.shape[0]):
        #     index_tuple_list = np.zeros(len(self.indices))
        #     index_tuple_list[0] = i
        #     index_tuple = tuple(index_tuple_list)
        #     print(index_tuple)
        #     vectors.append(components_tmp[i]*components_tmp.flatten()[i])
        #Maybe simply not a good idea.

def two_object_contraction(lorentz_object1,lorentz_object2):
    pass


#Dont have a good idea on this yet
#Takes a list of objects with Lorentz indices as arguments. The location
def general_contraction(lorentz_objects):
    pass

#2 -> 1 4vector add at vertices. The new particle is assumed to be on shell.
#Mass here is of little significance. Can override while knowing the real particle mass.
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

