import Lorentz_Transformation as lt
import LorentzGenerators as lg


class Indexed_L4v(lt.Lorentz4vector):
    def __init__(self,components,name='nothing',mass=None,index_name = 'mu'):
        lt.Lorentz4vector.__init__(self,components=components,mass=mass,name=name)
        self.index_name = index_name





