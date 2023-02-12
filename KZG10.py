from typing import List, NamedTuple, Tuple, Union
import json
from functools import reduce
from py_ecc.fields import bn128_FQ as FQ
from py_ecc import bn128 as curve
from py_ecc.typing import Field
from random import randint

G1_POINTS = []

class KZG10(object):

    def __init__(self) -> None:
        with open('taug1_65536.json', 'r') as file:
            g1_points = json.load(file)
        file.close()
        
        for point in g1_points:
            G1_POINTS.append(format_to_FQ(int(point[0], base=16), int(point[1], base=16)))

        self.F = GF(curve.curve_order)
        
    def evalPolyAt(self, coefficients: List[Field], index: Field):
        result = coefficients[0]
        for c_i in coefficients[1:]:
            result += c_i * index
            index = index*index
        return result
    
    def generate_coeffs(self, amount: int):
        return [self.F.random() for _ in range(amount-1)]
    
    def generate_proof(self, coefficients: List[Field], index: Field):
        n = len(coefficients)
        F = GF(curve.curve_order)
        y_powers = [F(1)]
        
        for i in range(n):
            y_powers.append(y_powers[-1] * index)

        result = None
        for i in range(0, n-1):
            for j in range(i, -1, -1):
                a = G1_POINTS[i]
                b = y_powers[i-j]
                c = coefficients[i+1]
                term = curve.multiply(a, int(b*c))
                result = term if result is None else curve.add(result, term)
        return result
    
    def get_index_x(self, index: int):
        return self.F(index)
    

    """
	Copute commitment to the evaluation of the polynomial given the coefficients
	"""
    def generate_commitment(self, coefficients: List[Field]):
        return reduce(curve.add, [curve.multiply(G1_POINTS[i], int(c_i))
							        for i, c_i in enumerate(coefficients)])

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
    for i in range(3):
        print(G1_POINTS[i])

if __name__ == "__main__":
    main()