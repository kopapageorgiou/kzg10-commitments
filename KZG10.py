from typing import List, NamedTuple, Tuple, Union
from math import ceil, log2
from random import randint
from functools import reduce
import operator
from py_ecc import bn128 as curve
from BlockchainLib import *
from BLS import *

"""
Implementation of PolyCommit_{DL} from:
"Constant-Size Commitments to Polynomials and Their Applications"
 - Aniket Kate, Gregory M. Zaverucha, and Ian Goldberg
 - Max Planck Institute for Software Systems (MPI-SWS)
 - https://www.cypherpunks.ca/~iang/pubs/PolyCommit-AsiaCrypt.pdf
 - http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf (extended version)
Section 3.2
"""



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

PointG1 = Tuple[curve.FQ, curve.FQ]
PointG2 = Tuple[curve.FQ2, curve.FQ2]


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


class TrustedSetup(NamedTuple):
	F: GF
	t: int
	g1_powers: List[PointG1]
	g2_powers: List[PointG2]
	alpha_powers: List[Field]

	@classmethod
	def generate(cls, F: GF, t: int, g1andg2: bool):
		"""
		Simulate the trusted setup
		Generate a random `alpha`
		Then return `t` powers of `alpha` in G1, and the representation of `alpha` in G2
		 g, g^a, ... g^{a^t}
		"""
		alpha = F.random()
		alpha_powers = [F(1)]
		g1_powers = [curve.G1]
		g2_powers = [curve.G2]
		for i in range(t+1):
			alpha_powers.append(alpha_powers[-1] * alpha)
			if g1andg2:
				g1_powers.append(curve.multiply(g1_powers[-1], int(alpha)))
				g2_powers.append(curve.multiply(g2_powers[-1], int(alpha)))
		
		return cls(F, t, g1_powers, g2_powers, alpha_powers)


def polynomial(x: Field, coeffs: List[Field]):
	result = coeffs[0]
	for c_i in coeffs[1:]:
		result += c_i * x
		x = x*x
	return result


def CommitProduct(PK: TrustedSetup, coeff: List[Field]):
	"""
	XXX: unsure if we need this, but it looks useful
	Computes commitment 'C' at `a`, where `a` is part of the trusted setup
		C = reduce(operator.mul, [(a**j) * c_j for j, c_j in enumerate(coeff)])
	For example:
		sage: factor((a^0 * phi0) * ((a^1) * phi1) * ((a^2) * phi2) * ((a^3)*phi3) * ((a^4)*phi4))
		a^10*(phi0*phi1*phi2*phi3*phi4)
	Where the element 'k' is the n-th triangle number
	"""
	n = len(coeff)
	k = (n*(n+1))//2  # equivalent to: reduce(operator.add, range(n))
	element = PK.g1_powers[k]
	product = PK.F(1)
	for c_i in coeffs:
		product *= c_i
	return curve.multiply(element, product)


def CommitSumTrusted(PK: TrustedSetup, coeff: List[Field]):
	return reduce(operator.add, [PK.alpha_powers[i] * c_i for i, c_i in enumerate(coeff)])


def CommitSum(PK: TrustedSetup, coeff: List[Field]):
	"""
	Copute commitment to the evaluation of a polynomial with coefficients
	At `x`, where `x` is part of the trusted setup
	"""
	return reduce(curve.add, [curve.multiply(PK.g1_powers[i], int(c_i))
							  for i, c_i in enumerate(coeff)])


def CommitRemainder(PK: TrustedSetup, y: Field, coeff: List[Field]):
	# TODO: implement
	"""
	f(x) = phi(x)
sage: f
x |--> c2*x^2 + c1*x + c0
sage: f.maxima_methods().divide((x-a) * (x-b))[1]
-a*b*c2 + (a*c2 + b*c2 + c1)*x + c0
sage: f.maxima_methods().divide((x-a))[1]
a^2*c2 + a*c1 + c0
sage: f.maxima_methods().divide((x-a) * (x-b) - (x-c))[1]
-a*b*c2 - c*c2 + (a*c2 + b*c2 + c1 + c2)*x + c0
	"""
	pass


def CommitDivisionTrusted(PK: TrustedSetup, y: Field, coeff: List[Field]):
	n = len(coeff)
	y_powers = [PK.F(1)]
	for i in range(n):
		y_powers.append(y_powers[-1] * y)

	result = PK.F(0)
	for i in range(0, n-1):
		for j in range(i, -1, -1):
			a = PK.alpha_powers[j]
			b = y_powers[i-j]
			c = coeff[i+1]
			result += a*b*c
	return result
	"""
	return reduce(operator.add, [PK.alpha_powers[j] * (y_powers[i-j] * coeff[i+1])
							     for j in range(i, -1, -1)
							     for i in range(0, n-1)])
	"""


def CommitDivision(PK: TrustedSetup, y: Field, coeff: List[Field]):
	"""
	Compute commitment to the division: `(phi(x) - phi(y)) / (x - y)`
	Using the trusted setup secret `x`
	Example:
		sage: poly4 = lambda k: c0 + ((k^1)*c1) + ((k^2)*c2) + ((k^3)*c3)
		sage: factor((poly4(x) - poly4(y)) / (x - y))
		c3*x^2 + c3*x*y + c3*y^2 + c2*x + c2*y + c1
	TODO: number of multiplications can be reduced significantly
	      number of additions can also be reduced significantly
	"""
	n = len(coeff)
	y_powers = [PK.F(1)]
	for i in range(n):
		y_powers.append(y_powers[-1] * y)

	result = None
	for i in range(0, n-1):
		for j in range(i, -1, -1):
			a = PK.g1_powers[j]
			b = y_powers[i-j]
			c = coeff[i+1]
			term = curve.multiply(a, int(b*c))
			result = term if result is None else curve.add(result, term)
	return result


def Prove():
	F = GF(curve.curve_order)
	coeff = [F.random() for _ in range(3)]
	PK = TrustedSetup.generate(F, len(coeff), True)

	"""
	The pairing equation works to verify that the witness and evaluation
	match the committed polynomial.
	sage:
		var('x y c0 c1 c2 c3')
		phi = lambda k: c0 + ((k^1)*c1) + ((k^2)*c2)
		psi = lambda a, b: ((phi(a) - phi(b)) / (a - b))
		psi(x, y) * (x-y) == phi(x) - phi(y)
		(psi3(x, y) * (x-y)) + phi(y) == phi(x)
	sage: psi(x, y) * (x-y)
	c2*x^2 - c2*y^2 + c1*x - c1*y
	sage: psi(x, y)
	(c2*x^2 - c2*y^2 + c1*x - c1*y)/(x - y)
	sage: factor(psi(x, y))
	c2*x + c2*y + c1
	sage: phi(x) - phi(y)
	c2*x^2 - c2*y^2 + c1*x - c1*y
	sage: factor(phi(x) - phi(y))
	(c2*x + c2*y + c1)*(x - y)
	"""

	# Verify with trusted information
	x = PK.alpha_powers[1]
	phi_at_x = polynomial(x, coeff)
	assert phi_at_x == CommitSumTrusted(PK, coeff)

	i = F(3)
	phi_at_i = polynomial(i, coeff)
	a = polynomial(x, coeff) - phi_at_i
	b = a / (x - i)

	psi_i_at_x = CommitDivisionTrusted(PK, i, coeff)
	assert psi_i_at_x == b
	assert psi_i_at_x * (x-i) == phi_at_x - phi_at_i


	# Then make commitment without access to trusted setup secrets
	g1_phi_at_x = CommitSum(PK, coeff)  # Commit to polynomial

	# Commit to an evaluation of the same polynomial at i
	i = F.random()  # randomly sampled i
	phi_at_i = polynomial(i, coeff)	
	g1_psi_i_at_x = CommitDivision(PK, i, coeff)

	# Compute `x - i` in G2
	g2_i = curve.multiply(curve.G2, int(i))
	g2_x_sub_i = curve.add(PK.g2_powers[1], curve.neg(g2_i)) # x-i

	# Verifier
	g1_phi_at_i = curve.multiply(curve.G1, int(phi_at_i))
	g1_phi_at_x_sub_i = curve.add(g1_phi_at_x, curve.neg(g1_phi_at_i))
	"""
	devPapas
	"""
	contract = smartContract()
	tx = contract.verify(formatG1(g1_phi_at_x_sub_i),formatG1(g1_phi_at_i), int(i), int(x))
	print(tx)
	#a = curve.pairing(g2_x_sub_i, g1_psi_i_at_x)
	#b = curve.pairing(curve.G2, curve.neg(g1_phi_at_x_sub_i))
	#ab = a*b
	#print(x)
	#print('ab', ab, ab == curve.FQ12.one())

if __name__ == "__main__":
	Prove()



def derp(n):
	"""
	Example:
	sage: poly5 = lambda k: c0 + ((k^1)*c1) + ((k^2)*c2) + ((k^3)*c3) + ((k^4)*c4) + ((k^5)*c5)
	sage: poly6 = lambda k: c0 + ((k^1)*c1) + ((k^2)*c2) + ((k^3)*c3) + ((k^4)*c4) + ((k^5)*c5)
	sage: poly5 = lambda k: c0 + ((k^1)*c1) + ((k^2)*c2) + ((k^3)*c3) + ((k^4)*c4)
	sage: poly4 = lambda k: c0 + ((k^1)*c1) + ((k^2)*c2) + ((k^3)*c3)
	sage: poly3 = lambda k: c0 + ((k^1)*c1) + ((k^2)*c2)
	sage: poly2 = lambda k: c0 + ((k^1)*c1)
	sage: poly2(x) - poly2(y)
	c1*x - c1*y
	sage: factor(poly2(x) - poly2(y))
	c1*(x - y)
	sage: factor((poly2(x) - poly2(y)) / (x - y))
	c1
	sage: factor((poly3(x) - poly3(y)) / (x - y))
	c2*x + c2*y + c1
	sage: factor((poly4(x) - poly4(y)) / (x - y))
	c3*x^2 + c3*x*y + c3*y^2 + c2*x + c2*y + c1
	sage: factor((poly5(x) - poly5(y)) / (x - y))
	c4*x^3 + c4*x^2*y + c4*x*y^2 + c4*y^3 + c3*x^2 + c3*x*y + c3*y^2 + c2*x + c2*y + c1
	sage: factor((poly6(x) - poly6(y)) / (x - y))
	c5*x^4 + c5*x^3*y + c5*x^2*y^2 + c5*x*y^3 + c5*y^4 + c4*x^3 + c4*x^2*y + c4*x*y^2 + c4*y^3 + c3*x^2 + c3*x*y + c3*y^2 + c2*x + c2*y + c1
	Emit this sequence:
	>>> derp(4)
	c[1] * x^0 * y^0
	c[2] * x^1 * y^0
	c[2] * x^0 * y^1
	c[3] * x^2 * y^0
	c[3] * x^1 * y^1
	c[3] * x^0 * y^2
	>>> derp(5)
	c[1] * x^0 * y^0
	c[2] * x^1 * y^0
	c[2] * x^0 * y^1
	c[3] * x^2 * y^0
	c[3] * x^1 * y^1
	c[3] * x^0 * y^2
	c[4] * x^3 * y^0
	c[4] * x^2 * y^1
	c[4] * x^1 * y^2
	c[4] * x^0 * y^3
	"""
	for i in range(0, n-1):
		for j in range(i, -1, -1):
			print(f'c[{i+1}] * x^{j} * y^{i-j}')
