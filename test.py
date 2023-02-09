#from py_eth_pairing import curve_add, curve_mul, pairing2, curve_negate
#from py_ecc.bn128 import bn128_field_elements
#from py_ecc.bn128 import bn128_curve
"""import sys
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
print(result)"""

#from py_ecc.bn128 import FQ, G1
from py_ecc.bn128 import bn128_curve as bn128
import sys
sys.path.append('./py_pairing-master/py_ecc/bn128')
from bn128_field_elements import FQ
from bn128_curve import add, G1
from bn128_curve import multiply

def polynomial_evaluation(coeffs, x):
    result = 0
    for i, coeff in enumerate(coeffs):
        result += coeff * (x**i)
    return result

def generate_commitment(coeffs, x):
    # Evaluate the polynomial
    result = polynomial_evaluation(coeffs, x)
    print(result)
    # Generate the commitment as a G1 point
    commitment = multiply(G1, result)

    return commitment

def generate_proof(coeffs, x, commitment):
    # Generate the proof using the given commitment and coefficients
    # (Note: this is a simplified example and not a complete implementation of a zero-knowledge proof)

    # Calculate the polynomial evaluation at x
    evaluation = polynomial_evaluation(coeffs, x)

    # Check that the evaluation and the commitment match
    if evaluation == commitment:
        # The proof is the commitment itself
        proof = commitment
    else:
        # The proof is None if the evaluation and commitment do not match
        proof = None

    return proof

# Example usage:
coeffs = [7261640105655323564317584541854666870196307670309556897762642496909946713003,
          5251439980146187779405380323454706999635699865129324559228435078229340579872,
          19093677638861175721480225560573248980575889224335987110361063931863808451363]
value = 1
#poly = polynomial_evaluation(coeffs, value)
#print(poly)
commitment = generate_commitment(coeffs, value)
print(commitment)
proof = generate_proof(coeffs, value, commitment)
print(proof)


