#from finitefield import GF
from numpy import array
from typing import List, NamedTuple, Tuple, Union
import json
from functools import reduce
from py_ecc.fields import bn128_FQ as FQ
from py_ecc import bn128 as curve
from py_ecc.typing import Field
from random import randint

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

    """def interp_poly(self, X, Y):
        poly = [[]]
        for j, y in enumerate(Y):
            Xe = X[:j] + X[j+1:]
            numer = reduce(lambda p, q: self.F(FQ(p))* self.F(FQ(q)), ([[1], K.sub([], x)] for x in Xe))
            denom = reduce(lambda x, y: K.mul(x, y), (K.sub(X[j], x) for x in Xe))
            poly = R.add(poly, R.mul(numer, [K.mul(y, K.inv(denom))]))
        return poly"""
    
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
    
# Define the finite field (in this case, GF(2^8))
F = GF(curve.curve_order)

# Define the data points
data_points = [(F(1), F(3)), (F(2), F(7)), (F(3), F(11)), (F(4), F(15)), (F(5), F(19))]

# Function to compute the Lagrange basis polynomial
def lagrange_basis_polynomial(x, i, data_points):
    numerator = F(1)
    denominator = F(1)
    for j, (xj, yj) in enumerate(data_points):
        if j == i:
            continue
        numerator *= x - xj
        denominator *= data_points[i][0] - xj
    return numerator / denominator

# Function to compute the Lagrange polynomial
def lagrange_polynomial(x, data_points):
    y = F(0)
    for i, (xi, yi) in enumerate(data_points):
        y += yi * lagrange_basis_polynomial(x, i, data_points)
    return y

# Compute the coefficients of the polynomial
x_values = [F(x) for x in range(6)]
A = array([[x**i for i in range(len(data_points))] for x in x_values])
b = array([lagrange_polynomial(x, data_points) for x in x_values])
print(b)
coefficients = array(F.inverse(A).dot(b), dtype=int)

print(coefficients)
