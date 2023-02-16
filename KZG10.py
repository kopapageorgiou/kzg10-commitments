from typing import List, Tuple
import json
from functools import reduce
from py_ecc.fields import bn128_FQ as FQ
from py_ecc.fields import bn128_FQ2 as FQ2
from py_ecc.bn128 import bn128_pairing as curvePairing
from py_ecc import bn128 as curve
from py_ecc.typing import Field
from random import randint
import numpy as np

G1_POINTS = []

SRS_G2_1 = (FQ2([
                0x04c5e74c85a87f008a2feb4b5c8a1e7f9ba9d8eb40eb02e70139c89fb1c505a9,
                0x21a808dad5c50720fb7294745cf4c87812ce0ea76baa7df4e922615d1388f25a
            ]),
            FQ2([
                0x2d58022915fc6bc90e036e858fbc98055084ac7aff98ccceb0e3fde64bc1a084,
                0x204b66d8e1fadc307c35187a6b813be0b46ba1cd720cd1c4ee5f68d13036b4ba
                
            ]))

class KZG10(object):

    def __init__(self, field = curve.curve_order) -> None:
        """
        Load G1 Points from tau json file
        and create GF with field= bn128.curve_order a.k.a BABYJUB_P (the same field with the one in solidity)
        """
        with open('taug1_65536.json', 'r') as file:
            g1_points = json.load(file)
        file.close()
        
        for point in g1_points:
            G1_POINTS.append(format_to_FQ(int(point[0], base=16), int(point[1], base=16)))
        self.field = field
        self.F = GF(field)
        
    def evalPolyAt(self, coefficients: List[Field], index: Field):
        """result = self.F(0)
        power_of_x = self.F(1)

        for coeff in coefficients:
            result = (result + (power_of_x * coeff))
            power_of_x = (power_of_x * index)

        return result"""
        result = self.F(0)
        n = len(coefficients)
        for i in range(n):
            result += coefficients[i] * (index ** (n-i-1))

        return result
    
    def generate_coeffs_random(self, amount: int):
        return [self.F.random() for _ in range(amount-1)]
    
    def generate_coeffs_for(self, values: List[int]):
        x_values = []
        val = []
        for x, y in enumerate(values):
            x_values.append(x)
            val.append(y)

        pol = self._lagrange_inter([self.F(x) for x in x_values], [self.F(y) for y in val])
        return(pol)
    
    def generate_proof(self, coefficients: List[Field], index: Field):
        quotientCoefficients = self._genQuotientPolynomial(coefficients, index)
        print("here", quotientCoefficients)
        return self.generate_commitment(quotientCoefficients)

    def _genQuotientPolynomial(self, coefficients: List[Field], xVal = Field):
        yVal = self.evalPolyAt(coefficients, xVal)
        x = [self.F(0), self.F(1)]
        res = self._divPoly(self._subPoly(coefficients, [yVal]), self._subPoly(x, [xVal]))[0]
        return res
    
    def get_index_x(self, index: int):
        return self.F(index)
    

    """
	Copute commitment to the evaluation of the polynomial given the coefficients
	"""
    def generate_commitment(self, coefficients: List[Field]):
        return reduce(curve.add, [curve.multiply(G1_POINTS[i], int(c_i))
							        for i, c_i in enumerate(coefficients)])

    def _subPoly(self, p1: List[Field], p2: List[Field]):
        degree = max(len(p1), len(p2))
        result = [self.F(0)] * degree
        for i in range(degree):
            if i < len(p1):
                result[i] += p1[i]
            if i < len(p2):
                result[i] -= p2[i]
        return result
    
    def _addPoly(self, p1, p2):
        if len(p1) < len(p2):
            p1 += [self.F(0)] * (len(p2) - len(p1))
        else:
            p2 += [self.F(0)] * (len(p1) - len(p2))
        
        # Perform element-wise addition of the coefficients
        result = [a + b for a, b in zip(p1, p2)]
        
        return result

    def _divPoly(self, numerator, denominator):
        degree_numerator = len(numerator) - 1
        degree_denominator = len(denominator) - 1
        if degree_numerator < degree_denominator:
            return [self.F(0)], numerator
        quotient = [self.F(0)] * (degree_numerator - degree_denominator + 1)
        remainder = numerator.copy()
        for i in range(degree_numerator - degree_denominator, -1, -1):
            quotient[i] = remainder[-1] / denominator[-1]
            for j in range(degree_denominator, -1, -1):
                remainder[i + j] -= quotient[i] * denominator[j]
            del remainder[-1]
        return quotient, remainder
    
    def _mulPoly(self, poly1, poly2):
    # Initialize an output polynomial with all coefficients set to zero
        output_poly = [self.F(0)] * (len(poly1) + len(poly2) - 1)
        # Multiply each term of poly1 with each term of poly2, and add the results to the output polynomial
        for i, coeff1 in enumerate(poly1):
            for j, coeff2 in enumerate(poly2):
                #print(coeff1, coeff2)
                output_poly[i+j] += coeff1 * coeff2

        return output_poly


    def _lagrange_interpolation(self, x_vals, y_vals):
        n = len(x_vals)
        coefficients = [self.F(0) for i in range(n)]
        for i in range(n):
            numerator = self.F(1)
            denominator = self.F(1)
            x_i = x_vals[i]
            y_i = y_vals[i]
            for j in range(n):
                if i == j:
                    continue
                numerator *= x_vals[j] - x_i
                denominator *= x_vals[j] - y_i
            coefficients[i] = y_i * numerator / denominator
        return coefficients
    
    def _lagrange_inter(self, x_vals, y_vals):
        L = [self.F(0)]
        k = len(x_vals)
        for j in range(k):
            lj = [y_vals[j]]
            for m in range(k):
                if m == j:
                    continue
                ljm = self._mulPoly([1, -x_vals[m]], [self._RECIPROCAL(int(x_vals[j] - x_vals[m]))])
                lj = self._mulPoly(lj, ljm)
            L = self._addPoly(L, lj)
        return L


    def _RECIPROCAL(self, a):
        if a == 0:
            return 0
        lm, hm = 1, 0
        low, high = a % self.field, self.field
        while low > 1:
            r = high//low
            nm, new = hm-lm*r, high-low*r
            lm, low, hm, high = nm, new, lm, low
        return lm % self.field

    def verify(self, commitment, proof, index, value):
        commitmentMinusA = curve.add(commitment, curve.neg(curve.multiply(curve.G1, int(value))))
        negProof = curve.neg(proof)
        indexMulProof = curve.multiply(proof, int(index))
        #return [commitmentMinusA, negProof, indexMulProof]

        a = curvePairing.pairing(curve.G2, curve.add(indexMulProof, commitmentMinusA))
        b = curvePairing.pairing(SRS_G2_1, negProof)
        ab = a*b
        print(ab == curve.FQ12.one())


class Field(object):
    def __init__(self, value, modulus: int):
        if isinstance(value, Field):
            value = value.v
        else:
            value = value % modulus
        self.v = value
        self.m = modulus

    def __eq__(self, other):
        if isinstance(other, int):
            return self.v == other % self.m
        elif isinstance(other, Field):
            return self.v == other.v and self.m == other.m
        else:
            raise ValueError(f'Cannot compare {self} with {other}')

    def __int__(self):
        return self.v

    def __add__(self, other):
        if isinstance(other, Field):
            other = other.v
        return Field(self.v + other, self.m)

    def __neg__(self):
        return Field(-self.v, self.m)

    def __mul__(self, other):
        if isinstance(other, Field):
            other = other.v
        return Field(self.v * other, self.m)

    def __repr__(self):
        return f"Field<{self.v}>"

    def __sub__(self, other):
        if isinstance(other, Field):
            other = other.v
        return Field(self.v - other, self.m)

    def __truediv__(self, other):
        return self * other.inverse()

    def __pow__(self, other):
        assert isinstance(other, int)
        return Field(pow(self.v, other, self.m), self.m)

    def inverse(self):
        return Field(pow(self.v, self.m-2, self.m), self.m)
   
class GF(object):
    def __init__(self, modulus: int):
        self.m = modulus

    def primitive_root(self, n: int):
        """
        Find a primitive n-th root of unity
            - http://www.csd.uwo.ca/~moreno//AM583/Lectures/Newton2Hensel.html/node9.html
            - https://crypto.stackexchange.com/questions/63614/finding-the-n-th-root-of-unity-in-a-finite-field
        """
        assert n >= 2
        x = 2
        while True:
            # x != 0, g = x^(q-1/n)
            # if g^(n/2) != 1, then it is a primitive n-th root
            g = pow(int(x), (self.m-1)//n, self.m)
            if pow(g, n//2, self.m) != 1:
                return self(g)
            x += 1

    def random(self):
        return self(randint(0, self.m-1))

    def __call__(self, value) -> Field:
        return Field(value, self.m)


    
def format_to_FQ(x_point: int, y_point: int):
    return (FQ(x_point), FQ(y_point))

def format_FQ_G1Point(data: Tuple[Field, Field]):
    x, y = data
    return (int(str(x)), int(str(y)))

def format_field_to_int(value: Field | List[Field]):
    if (type(value) == Field):
        return int(value)
    else:
        return [int(val) for val in value]
    
    
    

def main():
    kzg = KZG10(127)
    coeffs = kzg.generate_coeffs_for([5, 25, 125])
    commit = kzg.generate_commitment(coeffs)
    print(coeffs)
    x = kzg.get_index_x(1)
    y = kzg.evalPolyAt(coeffs, x)

    proof = kzg.generate_proof(coeffs, x)
    kzg.verify(commit, proof, x, y)

if __name__ == "__main__":
    main()
