import matplotlib
import numpy as np

def sampleExp(Lambda =2, maxCnt = 50000):
    ys =[]
    standardXaxis = []
    standardexp = []
    for i in range(maxCnt):
        u = np.random.random()
        y = -1/Lambda*np.log(1-u)
        ys.append(y)
    for i in range(1000):
        t = Lambda * np.exp(
            \]
        )