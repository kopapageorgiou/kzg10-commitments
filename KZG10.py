from typing import List, Tuple
import json
from functools import reduce
from py_ecc.fields import bn128_FQ as FQ
from py_ecc import bn128 as curve
from py_ecc.typing import Field
from random import randint

G1_POINTS = []

class KZG10(object):

    def __init__(self) -> None:
        """
        Load G1 Points from tau json file
        and create GF with field= bn128.curve_order a.k.a BABYJUB_P (the same field with the one in solidity)
        """
        with open('taug1_65536.json', 'r') as file:
            g1_points = json.load(file)
        file.close()
        
        for point in g1_points:
            G1_POINTS.append(format_to_FQ(int(point[0], base=16), int(point[1], base=16)))

        self.F = GF(curve.curve_order)
        
    def evalPolyAt(self, coefficients: List[Field], index: Field):
        result = self.F(0)
        power_of_x = self.F(1)

        for coeff in coefficients:
            result = (result + (power_of_x * coeff))
            power_of_x = (power_of_x * index)

        return result
    
    def generate_coeffs_random(self, amount: int):
        return [self.F.random() for _ in range(amount-1)]
    
    def generate_coeffs_for(self, values: List[int]):
        x_values = []
        val = []
        for x, y in enumerate(values):
            x_values.append(y)
            val.append(x)

        pol = self._lagrange_interpolation([self.F(x) for x in x_values], [self.F(y) for y in val])
        return(pol)
    
    def generate_proof(self, coefficients: List[Field], index: Field):
        quotientCoefficients = self._genQuotientPolynomial(coefficients, index)
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
    kzg = KZG10()
    coeffs = kzg.generate_coeffs_2([5, 25, 30, 80])
    print(coeffs)

if __name__ == "__main__":
    main()
