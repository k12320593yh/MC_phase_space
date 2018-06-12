import re as re

test_string = 'abc123456abcdefghelloworld[aaa]helloworld'
test_string2 = '<'+test_string
test_pattern = re.compile(r'<a')
test_pattern2 = re.compile(r'abc(?P<word>.*)abc')
print(re.sub(pattern=test_pattern,string=test_string,repl='!'))

# abc123456abcdefghelloworld[aaa]helloworld

print(re.sub(pattern=test_pattern,string=test_string2,repl='!'))

# !bc123456abcdefghelloworld[aaa]helloworld
# Pointy parenthesis has no further meaning in re.
# Which is great!

# test_string3 = 'e^6*(Pair[Momentum[k1],Momentum[k2]])'
# test_string3 = r'(e^6*((-16*Pair[Momentum[k1], Momentum[p2]]*Pair[Momentum[k2], Momentum[p1]]*Pair[Momentum[k3], Momentum[k3]])/(Pair[Momentum[k3], Momentum[k3]] - 2*Pair[Momentum[k3], Momentum[p1]] + Pair[Momentum[p1], Momentum[p1]])^2 + (64*sw2*Pair[Momentum[k1], Momentum[p2]]*Pair[Momentum[k2], Momentum[p1]]*Pair[Momentum[k3], Momentum[k3]])/(Pair[Momentum[k3], Momentum[k3]] - 2*Pair[Momentum[k3], Momentum[p1]] + Pair[Momentum[p1], Momentum[p1]])^2 - (128*sw2^2*Pair[Momentum[k1], Momentum[p2]]*Pair[Momentum[k2], Momentum[p1]]*Pair[Momentum[k3], Momentum[k3]])/(Pair[Momentum[k3], Momentum[k3]] - 2*Pair[Momentum[k3], Momentum[p1]] + Pair[Momentum[p1], Momentum[p1]])^2 - (16*Pair[Momentum[k1], Momentum[p1]]*Pair[Momentum[k2], Momentum[p2]]*Pair[Momentum[k3], Momentum[k3]])/(Pair[Momentum[k3], Momentum[k3]] - 2*Pair[Momentum[k3], Momentum[p1]] + Pair[Momentum[p1], Momentum[p1]])^2 + (64*sw2*Pair[Momentum[k1], Momentum[p1]]*Pair[Momentum[k2], Momentum[p2]]*Pair[Momentum[k3], Momentum[k3]])/(Pair[Momentum[k3], Momentum[k3]] - 2*Pair[Momentum[k3], Momentum[p1]] + Pair[Momentum[p1], Momentum[p1]])^2 - (128*sw2^2*Pair[Momentum[k1], Momentum[p1]]*Pair[Momentum[k2], Momentum[p2]]*Pair[Momentum[k3], Momentum[k3]])/(Pair[Momentum[k3], Momentum[k3]] - 2*Pair[Momentum[k3], Momentum[p1]] + Pair[Momentum[p1], Momentum[p1]])^2 + (32*Pair[Momentum[k1], Momentum[p2]]*Pair[Momentum[k2], Momentum[k3]]*Pair[Momentum[k3], Momentum[p1]])/(Pair[Momentum[k3], Momentum[k3]] - 2*Pair[Momentum[k3], Momentum[p1]] + Pair[Momentum[p1], Momentum[p1]])^2'
with open("m2.txt","r") as file_object:
    test_string3 = file_object.read()
    # print(st)

test_patternP = re.compile(r'Pair')
test_patternM =re.compile(r'Momentum')
test_patternb = re.compile(r'\[\[')
test_patternbb = re.compile(r'\]\]')
test_patternbbb = re.compile(r'\^')
test_patternbbbb = re.compile(r'\)\)\)')
test_pattern5 = re.compile('], \[')
test_string3 = re.sub(string=test_string3,pattern=test_patternP,repl='lt.fcc')
print(test_string3)
test_string3 = re.sub(string=test_string3,pattern=test_patternM,repl='')
print(test_string3)
test_string3 = re.sub(string=test_string3,pattern=test_patternb,repl='(')
print(test_string3)
test_string3 = re.sub(string=test_string3,pattern=test_patternbb,repl=')')
# print(test_string3)
test_string3 = re.sub(string=test_string3,pattern=test_patternbbb,repl='**')
# print(test_string3)
# test_string3 = re.sub(string=test_string3,pattern=test_patternbbbb,repl=')')
print(test_string3)
test_string3 = re.sub(string=test_string3,pattern=test_pattern5,repl=',')
print(test_string3)
new_file_object = open('m2n','w')
new_file_object.write(test_string3)
new_file_object.close()