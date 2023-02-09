#from py_eth_pairing import curve_add, curve_mul, pairing2, curve_negate
#from py_ecc.bn128 import bn128_field_elements
#from py_ecc.bn128 import bn128_curve
import sys
sys.path.append('./py_pairing-master/py_ecc/bn128')
from bn128_field_elements import FQ
from bn128_curve import add
from bn128_curve import multiply
from polynomial import *
from Constants import *

coeffs = [7261640105655323564317584541854666870196307670309556897762642496909946713003,
          5251439980146187779405380323454706999635699865129324559228435078229340579872,
          19093677638861175721480225560573248980575889224335987110361063931863808451363]

result = (FQ(0), FQ(0))
for i,c in enumerate(coeffs):
    result = add(result, multiply((FQ(SRS_G1_X[i]), FQ(SRS_G1_Y[i])), c))
print(result)