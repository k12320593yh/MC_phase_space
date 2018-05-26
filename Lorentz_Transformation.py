import numpy as np
import LorentzGenerators as lg
import scipy.linalg as sp

#A general Lorentz 4 vector. Just what you learned in peskin.
class Lorentz4vector(object):
    def __init__(self,components,mass = None,name = 'nothing'):
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
        self.__index_case = True

    def get_4_vector(self):
        return np.array(self.__components)

    def get_mass(self):
        return self.__mass

    def get_on_shell_ness(self):
        return self.__on_shell

    def lorentz_transformation(self,ltf):
        # tmp = np.matmul(lg.minkowski_metric,ltf.get_matrix_form())
        self.__components = np.matmul(ltf.get_matrix_form(),self.__components)
        # self.__index_case = not self.__index_case

    def lorentz_transformation_off_matrix(self,matrix):
        self.__components = np.matmul(matrix,self.__components)

    def lorentz_rot_toz(self):
        rot_mat = find_toz_rotation(self.__components)
        # tmp = np.matmul(lg.minkowski_metric,rot_mat)
        self.__components = np.matmul(rot_mat,self.__components)
        # self.transformation = np.matmul(rot_mat,self.transformation)
        # self.__index_case = not self.__index_case

    def zrapidity(self):
        self.__zrapidity = find_boost_to_rest(self)

    def get_zrapidity(self):
        return self.__zrapidity

    #Must be done after rotate into z direction and
    def z_boost_to_rest_frame(self):
        boost_mat = boost_matrix(rapidities=np.array([0,0,-self.__zrapidity]))
        # print(boost_mat)
        # print(simple_z_boost(-self.__zrapidity))
        tmp = np.matmul(lg.minkowski_metric,boost_mat)
        self.__components = np.matmul(tmp,self.__components)
        self.transformation = np.matmul(boost_mat,self.transformation)
        self.__index_case = not self.__index_case

    # def four_vector_add(self,l4v):
    #     self.__components += l4v.get_4_vector()
    #     self.__mass = np.sqrt(self.__components[0]**2 - np.sum(self.__components[i]**2 for i in range(1,4)))
        #Propagator mass. Should be able to specific if the process is known.

    # def contraction(self,l4v):

#This might look nasty but actually as it seems.


#2 -> 1 4vector add at vertices. The new particle is assumed to be on shell.
#Mass here is of little significance. Can override while knowing the real particle mass.
def four_vector_add(l4v1,l4v2):
    summ = l4v1.get_4_vector()+l4v2.get_4_vector()
    #p^2 = m^2 for every on shell particle
    pmass = np.sqrt(summ[0]**2 - np.sum(summ[i]**2 for i in range(1,4)))
    return Lorentz4vector(components=summ,name='new vector',mass=pmass)


#4vector contracion. Returning a lorentz invariant scalar.
def four_vector_contraction(l4v1,l4v2):
    v1 = l4v1.get_4_vector()
    v2 = l4v2.get_4_vector()
    return v1[0]*v2[0] - np.sum(v1[i]*v2[i] for i in range(1,4))

def four_vector_contraction_component(l4v1,l4v2):
    return l4v1[0]*l4v2[0]-np.sum(l4v1[i]*l4v2[i] for i in range(1,4))

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
        if parameter != None:
            self.__parameter = parameter
            self.__matrix_form = general_matrix_form(self.__parameter)
        elif matrix != None:
            self.__matrix_form = matrix


    def show_matrix_form(self):
        print(self.__matrix_form)

    def get_matrix_form(self):
        return self.__matrix_form

    #backdoor measure to overwrite the matrix of a LorentzTransFormation object.
    #If you're going to use this method, initial the stance with a straight 0 parameter ndarray.
    def set_matrix_form(self,matrix):
        self.__matrix_form = matrix

