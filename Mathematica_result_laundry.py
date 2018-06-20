import numpy as np
import re

def wash(file):
    with open(file, "r") as file_object:
        st2 = file_object.read()

    pat1 = re.compile(r"Pair\[")
    st2 = re.sub(pattern=pat1, string=st2, repl=r'lt.fcc(')
    pat2 = re.compile(r"Momentum\[")
    st2 = re.sub(pattern=pat2, string=st2, repl=r'')
    pat3 = re.compile(r'], ')
    st2 = re.sub(pattern=pat3, string=st2, repl=',')
    pat4 = re.compile(r']]')
    st2 = re.sub(pattern=pat4, string=st2, repl=')')
    pat5 = re.compile('\^')
    st2 = re.sub(pattern=pat5, string=st2, repl='**')

    return st2