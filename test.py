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
import sys, os, pickle
sys.path.append('./py_pairing-master/py_ecc/bn128')
from bn128_field_elements import FQ
from bn128_curve import add, G1
from bn128_curve import multiply
import py_ecc.bn128 as curve
import numpy as np
from Constants import *
import random
from math import ceil, log2
import galois

BABYJUB_P = 21888242871839275222246405745257275088548364400416034343698204186575808495617

def as_bits_bytes(n):
	n_bits = ceil(log2(n))
	n_bytes = (n_bits + (8 - (n_bits % 8))) // 8
	return n_bits, n_bytes

FIELD_MODULUS = curve.field_modulus
FIELD_MODULUS_bits, FIELD_MODULUS_bytes = as_bits_bytes(FIELD_MODULUS)

CURVE_ORDER = curve.curve_order
CURVE_ORDER_bits, CURVE_ORDER_bytes = as_bits_bytes(CURVE_ORDER)

assert FIELD_MODULUS_bytes == 32
assert CURVE_ORDER_bytes == 32

def lagrange_interpolation(x_values, y_values):
    n = len(x_values)
    coeffs = []
    for i in range(n):
        x = x_values[i]
        y = y_values[i]
        l = 1
        for j in range(n):
            if j != i:
                l *= (x - x_values[j]) / (x_values[i] - x_values[j])
        coeffs.append(int(y * l))
    return np.array(coeffs)

def polynomial_interpolation():
    
    coeffs = lagrange_interpolation(SRS_G1_X, SRS_G1_Y)
    def polynomial(x):
        result = 0
        for i, coeff in enumerate(coeffs):
            result += coeff * (x**i)
        return result
    return polynomial

poly = polynomial_interpolation()

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

def generate_commitment_2(coeffs):
    print(coeffs)
    result = curve.add(curve.multiply((FQ(SRS_G1_X[0]), FQ(SRS_G1_Y[0])), coeffs[0]),
                        curve.multiply((FQ(SRS_G1_X[1]), FQ(SRS_G1_Y[1])), coeffs[1]))
    for i in range(2, len(coeffs)):
        temp = (FQ(SRS_G1_X[i]), FQ(SRS_G1_Y[i]))
        multemp = curve.multiply(temp, coeffs[i])
        result = curve.add(result, multemp)
        #print(SRS_G1_X[i], SRS_G1_Y[i])
        #result = curve.add(result, curve.multiply((FQ(SRS_G1_X[i]), FQ(SRS_G1_Y[i])), coeffs[i]))
    return result

def testAdd(coeffs):
    result = curve.add((FQ(SRS_G1_X[0]), FQ(SRS_G1_Y[0])), (FQ(SRS_G1_X[1]), FQ(SRS_G1_Y[1])))
    #for i in range(2, len(sigs)):
    #    aggregated = add(aggregated, sigs[i])
    #print(aggregated)
    #result = (FQ(0), FQ(0))
    #print(type(result[0]))
    for i in range(2, len(coeffs)):
        #print(SRS_G1_X[i], SRS_G1_Y[i])
        result = curve.add(result, (FQ(SRS_G1_X[i]), FQ(SRS_G1_Y[i])))
    return result

def testMul(coeffs):
    result = (FQ(SRS_0), FQ(SRS_0))
    for i in range(len(coeffs)):
        result = curve.multiply((FQ(SRS_G1_X[i]), FQ(SRS_G1_Y[i])), coeffs[i])
    return result
def generate_proof(coeffs, x, commitment):
    # Generate the proof using the given commitment and coefficients
    # (Note: this is a simplified example and not a complete implementation of a zero-knowledge proof)

    # Calculate the polynomial evaluation at x
    evaluation = polynomial_evaluation(coeffs, x)

    # Check that the evaluation and the commitment match

    return evaluation

def genProof(coeffs, x):
    print(type(x))
    print(coeffs)
    (new_vals, yval) = _genQuotientPolynomial(coeffs, x)
    print(new_vals)
    return (generate_commitment_2(new_vals), yval)

def _genQuotientPolynomial(coeffs, xval):
    if os.path.isfile("babyjub.field"):
        pass
    field = galois.Field(BABYJUB_P)
    polyn = galois.Poly(coeffs, field=field)
    print(polyn)
    yval = polyn(at=xval, field=field)
    print(yval)
    y = galois.Poly([yval], field=field)
    print(y)
    x = galois.Poly([0,1], field=field)
    print(x)
    z = galois.Poly([xval], field=field)
    print(z)
    result = galois.gcd(polyn-y,x-z)
    print("res=", result)
    print(result.coefficients())
    return([int(co) for co in result.coefficients()], int(yval)) 
def generate_random_x():
    random_integer = random.randint(0, 2**256 - 1)
    return FQ(random_integer)
# Example usage:
def main():
    coeffs = [2, 20572311769411586958468171033108896883420015526259878468433894867106118238149, 12142410488737504677019842007390042101436112662914794898876151082116442432673]
    value = 20874104960467819026323991510096759805217606586237124613779825101767616141755
    #poly = polynomial_evaluation(coeffs, value)
    #print(poly)
    commitment = generate_commitment(coeffs, value)
    print(commitment)
    proof = generate_proof(coeffs, value, commitment)
    print(proof)
    print(generate_commitment_2(coeffs))
    print(genProof(coeffs, value))

if __name__=="__main__":
    main()

