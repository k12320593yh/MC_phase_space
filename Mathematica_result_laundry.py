import numpy as np
import re

st = "-16*Pair[Momentum[k1], Momentum[p2]]*Pair[Momentum[k2], Momentum[p1]]"

def wash(file):
    with open(file, "r") as file_object:
        st2 = file_object.read()

    pat1 = re.compile(r"Pair")
    st2 = re.sub(pattern=pat1, string=st2, repl=r'lt.fcc')
    pat2 = re.compile(r"Momentum")
    st2 = re.sub(pattern=pat2, string=st2, repl=r'')
    pat3 = re.compile(r']')
    st2 = re.sub(pattern=pat3, string=st2, repl=')')
    pat4 = re.compile(r'\[')
    st2 = re.sub(pattern=pat4, string=st2, repl='(')
    pat5 = re.compile('\^')
    st2 = re.sub(pattern=pat5, string=st2, repl='**')

    return st2

print(wash('m2.txt'))
#This wash function passed the mathematica validation with a 1/1000 difference, which can be considered as 'good enough'.