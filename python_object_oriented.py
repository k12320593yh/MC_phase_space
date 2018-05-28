class Excretment(object):
    def __init__(self):
        self.a = 1

class Shit(Excretment):
    def __init__(self):
        Excretment.__init__(self)

a = Shit()
print(a.a)