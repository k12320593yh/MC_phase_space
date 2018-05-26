import numpy as np
import re

teststring = 'abcdefg123456abc'
teststring2 = "(2 (Overscript[k1, _]\[CenterDot]Overscript[p1, _]) (Overscript[k2, _]\[CenterDot]Overscript[p2, _]) Overscript[p2, _]^2 e^6)/(cw2^2 mz2^2 (Overscript[k3, _]^2-2 (Overscript[k3, _]\[CenterDot]Overscript[p2, _])+Overscript[p2, _]^2)^2)"
print(re.search(r'',teststring).group())
