import numpy as np
import re

with open("m2.txt","r") as file_object:
    st = file_object.read()
    print(st)
# teststring = 'abcdefg123456abc'
# teststring2 = "(2 (Overscript[k1, _]\[CenterDot]Overscript[p1, _]) (Overscript[k2, _]\[CenterDot]Overscript[p2, _]) Overscript[p2, _]^2 e^6)/(cw2^2 mz2^2 (Overscript[k3, _]^2-2 (Overscript[k3, _]\[CenterDot]Overscript[p2, _])+Overscript[p2, _]^2)^2)"
# print(re.search(r'',teststring))
# pattern = re.compile(r"smart")
# line = "I'm very smart"
# matchObj = re.search(pattern=pattern,string=line)
# # print(matchObj.group())
# line1 = "(e^6*((-16*Pair[Momentum[k1], Momentum[p2]]*Pair[Momentum[k2], Momentum[p1]]*Pair[Momentum[k3], Momentum[k3]])/(Pair[Momentum[k3], Momentum[k3]] - 2*Pair[Momentum[k3], Momentum[p1]] + Pair[Momentum[p1], Momentum[p1]])^2 + (64*sw2*Pair[Momentum[k1],"
# pattern = re.compile(r"Pair")
# pattern1 = re.compile(r"")
# print(re.search(pattern=pattern,string=line1).group())
# test_pattern = re.compile("Pair\\[Momentum\\[..\\], Momentum\\[..\\]\\]")
# print(re.sub(pattern=test_pattern,repl="shit",string=line1))
