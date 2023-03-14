from typing import List, NamedTuple, Tuple, Union
from random import randint
from math import ceil, log2
from functools import reduce
import operator
from py_ecc import bn128 as curve
from BlockchainLib import *
from BLS import *
import sys
from py_ecc.fields import bn128_FQ2 as FQ2
from py_ecc.typing import Point2D, Field
from KZG10 import *
import random
import time


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
		print("a=", alpha)
		#alpha = F(21888242871839275222246405745257275088548364400416034343698204186575808495617)
		alpha_powers = [F(1)]
		g1_powers = [curve.G1]
		g2_powers = [curve.G2]
		for i in range(t):
			#alpha_powers.append(alpha_powers[-1] * alpha)
			alpha_powers.append(alpha_powers[-1] * alpha)
			if g1andg2:
				g1_powers.append(curve.multiply(g1_powers[-1], int(alpha)))
				g2_powers.append(curve.multiply(g2_powers[-1], int(alpha)))
		print("g1_powers= ", g1_powers)
		
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
	print(product)
	for c_i in coeff:
		product *= c_i
	return curve.multiply(element, product)


def CommitSumTrusted(PK: TrustedSetup, coeff: List[Field]):
	return reduce(operator.add, [PK.alpha_powers[i] * c_i for i, c_i in enumerate(coeff)])


def CommitSum(coeff: List[Field]):
	"""
	Copute commitment to the evaluation of a polynomial with coefficients
	At `x`, where `x` is part of the trusted setup
	"""
	print("here2",reduce(curve.add, [curve.multiply((FQ(SRS_G1_X[i]), FQ(SRS_G1_Y[i])), int(c_i))
							  for i, c_i in enumerate(coeff)]))
	return reduce(curve.add, [curve.multiply((FQ(SRS_G1_X[i]), FQ(SRS_G1_Y[i])), int(c_i))
							  for i, c_i in enumerate(coeff)])

def mul_scalar(point, k: int):


	if k < 0:
	# k * point = -k * (-point)
		return mul_scalar(-k, point_neg(point))

	result = None
	addend = point

	while k:
		if k & 1:
		# Add.
			result = add(result, addend)

	# Double.
		addend = add(addend, addend)

		k >>= 1

		assert is_on_curve(result, FQ(3))
		#print(result)

	return result

def point_neg(point):
    """Returns -point."""
    assert is_on_curve(point)

    if point is None:
        # -0 = 0
        return None

    x, y = point
    result = (x, -y % curve.p)

    assert is_on_curve(result)

    return result

def point_add(point1, point2):
    """Returns the result of point1 + point2 according to the group law."""
    assert is_on_curve(point1, FQ(3))
    assert is_on_curve(point2, FQ(3))

    if point1 is None:
        # 0 + point2 = point2
        return point2
    if point2 is None:
        # point1 + 0 = point1
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        # point1 + (-point1) = 0
        return None

    if x1 == x2:
        # This is the case point1 == point2.
        m = (3 * x1 * x1 + curve.a) * inverse_mod(2 * y1, curve.p)
    else:
        # This is the case point1 != point2.
        m = (y1 - y2) * inverse_mod(x1 - x2, curve.p)

    x3 = m * m - x1 - x2
    y3 = y1 + m * (x3 - x1)
    result = (x3 % curve.p, -y3 % curve.p)

    assert is_on_curve(result)

    return result

def inverse_mod(k, p):
    """Returns the inverse of k modulo p.
    This function returns the only integer x such that (x * k) % p == 1.
    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError('division by zero')

    if k < 0:
        # k ** -1 = p - (-k) ** -1  (mod p)
        return p - inverse_mod(-k, p)

    # Extended Euclidean algorithm.
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, k

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    gcd, x, y = old_r, old_s, old_t

    assert gcd == 1
    assert (k * x) % p == 1

    return x % p


# Functions that work on curve points #########################################

def commit(coeffs):
	result = (FQ(0), FQ(0))
	#print("coeffs",coeffs)
	for i,c in enumerate(coeffs):
		
		#print(i,c)
		result = add(result, curve.multiply((FQ(int(SRS_G1_X[i],base=16)), FQ(int(SRS_G1_Y[i],base=16))), c))
		#print(FQ(int(SRS_G1_X[i],base=16)), FQ(int(SRS_G1_Y[i],base=16)))
	return result
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


def CommitDivision(y: Field, coeff: List[Field]):
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
	F = GF(curve.curve_order)
	n = len(coeff)
	y_powers = [F(1)]
	
	for i in range(n):
		y_powers.append(y_powers[-1] * y)
		print(y_powers)

	result = None
	for i in range(0, n-1):
		for j in range(i, -1, -1):
			a = (FQ(SRS_G1_X[i]), FQ(SRS_G1_Y[i]))
			b = y_powers[i-j]
			c = coeff[i+1]
			term = curve.multiply(a, int(b*c))
			result = term if result is None else curve.add(result, term)
	return result


def Prove():
	contract = smartContract()
	F = GF(curve.curve_order)
	number_of_coeffs = 3
	coeff = [F.random() for _ in range(number_of_coeffs-1)]
	int_coeff = [int(c) for c in coeff]
	#print(coeff)
	#PK = TrustedSetup.generate(F, len(coeff), True)

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
	"""
	# ? Excluded
	# Verify with trusted information
	x = PK.alpha_powers[1]
	phi_at_x = polynomial(x, coeff)
	#print("phi_at_x",phi_at_x)
	assert phi_at_x == CommitSumTrusted(PK, coeff)
	"""
	
	"""
	# ? Excluded
	i = F(3)
	phi_at_i = polynomial(i, coeff)
	#print("phi_at_i",phi_at_i)
	a = polynomial(x, coeff) - phi_at_i
	b = a / (x - i)

	psi_i_at_x = CommitDivisionTrusted(PK, i, coeff)
	#print("psi_i_at_x",psi_i_at_x)
	assert psi_i_at_x == b
	assert psi_i_at_x * (x-i) == phi_at_x - phi_at_i
	"""


	# Then make commitment without access to trusted setup secrets
	#g1_phi_at_x = CommitSum(coeff)  # Commit to polynomial
	#print("g1_phi_at_x", g1_phi_at_x)
	# Commit to an evaluation of the same polynomial at i
	index_x = F.random()  #! WE GOT i 
	print("index-x =", index_x)
	value_y = polynomial(index_x, coeff)	
	print("y value =", value_y) #! CORRECT
	commit = CommitSum(coeff)
	print("commitment =", commit) #! CORRECT
	tx = contract.commit(int_coeff)
	print("contract", tx)
	print(tx == (formatG1_FQ(commit)))
	proof = CommitDivision(index_x, coeff)
	print("proof =", proof)
	tx = contract.verify(formatG1_FQ(commit), formatG1_FQ(proof), int(index_x), int(value_y))
	print("contract", tx)
	#tx = contract.commit([int(c) for c in coeff])
	#print("contract", tx)

	# Compute `x - i` in G2
	#g2_i = curve.multiply(curve.G2, int(i))
	#print("g2_i", g2_i)
	#g2_x_sub_i = curve.add(PK.g2_powers[1], curve.neg(g2_i)) # x-i
	#print("g2_x_sub_i", g2_x_sub_i)

	# Verifier
	#g1_phi_at_i = curve.multiply(curve.G1, int(index_x))
	#print("g1_phi_at_i", g1_phi_at_i)
	#g1_phi_at_x_sub_i = curve.add(commit, curve.neg(g1_phi_at_i))
	#print("g1_phi_at_x_sub_i", g1_phi_at_x_sub_i)
	"""
	devPapas
	"""
	#coeffs = [int(c) for c in coeff]
	#print(coeffs)
	
	#print(formatG1(curve.neg(g1_phi_at_x_sub_i)),"\n",formatG2(curve.G2),"\n", formatG1(g1_psi_i_at_x),"\n", formatG2(g2_x_sub_i))
	#tx = contract.pairing(formatG1(curve.neg(g1_phi_at_x_sub_i)), formatG2(curve.G2), formatG1(g1_psi_i_at_x), formatG2(g2_x_sub_i))
	#print("pairing test =", tx)
	#tx = contract.evalPolyAt(coeffs, int(i))
	#print("contract", tx)
	#print("custom",g1_phi_at_x)
	#tx = contract.commit(coeffs)
	#print("contract", tx)
	
	#tx = contract.pairing(	formatG1(g1_psi_i_at_x),
	#						formatG2(curve.G2),
	#						formatG1(curve.neg(g1_phi_at_x_sub_i)),
	#	       				formatG2(g2_x_sub_i)
	#						)
	#print("contract", tx)
	#a = curve.pairing(g2_x_sub_i, g1_psi_i_at_x)
	#b = curve.pairing(curve.G2, curve.neg(g1_phi_at_x_sub_i))
	#ab = a*b
	#print(a)
	#print(b)
	#print('ab', ab, ab == curve.FQ12.one())


def genereateRandoms():
	unique_numbers = set()
	# generate random numbers until we have 100 unique ones
	while len(unique_numbers) < 350:
		unique_numbers.add(random.randint(1, 1000000))

	# convert the set to a list and sort it
	random_numbers = sorted(list(unique_numbers))
	return random_numbers


def testRun():
	#contract = smartContract()

	kzg = KZG10()
	
	start_time = time.time()  # get the current time in seconds

	coeffs = kzg.generate_coeffs_for([5, 66, 120, 212])

	end_time = time.time()  # get the current time again

	elapsed_time = end_time - start_time  # calculate the elapsed time in seconds

	print("Elapsed time:", elapsed_time, "seconds")

	start_time = time.time()  # get the current time in seconds

	x = np.arange(0, 200, 1)
	#y = random_numbers
	coefs2 = kzg.generate_coeffs_for2([5, 66, 120, 212])

	#a_s = divided_diff(x, random_numbers)[0, :]

	#x_new = np.arange(0, 200, 1)
	#y_new = newton_poly(a_s, x, x_new)
	
	elapsed_time = end_time - start_time  # calculate the elapsed time in seconds

	print("Elapsed time:", elapsed_time, "seconds")
	print("Coefs : ", coeffs)
	print("Coefs2 : ", coefs2[0])
	
	#coeffs = kzg.generate_coeffs(10)
	#print("coefficients = \n", x_new)
	#print("coefficients = \n", y_new)
	#index_x = kzg.get_index_x(1)  					#! Choosing index_x
	#print("\nIndex x:\n", index_x)
	#print("index-x =", index_x)
	#value_y = kzg.evalPolyAt(coeffs, index_x)		#! Evaluating polynomial at index_x to get y_value
	#print("\ny value:\n", value_y)



	#tx = contract.evalPolyAt2(format_field_to_int(coeffs), format_field_to_int(index_x))
	#print(tx)
	#print("\nIs the y value equal with the value from the chain?", format_field_to_int(value_y) == tx)
	#commit = kzg.generate_commitment(coeffs)		#! Generating commitment
	#print("\nCommitment: \n", commit)
	#tx = contract.commit(format_field_to_int(coeffs))
	#print("\nIs the commitment equal with the commitment from the chain?", tx == format_FQ_G1Point(commit))
	#proof = kzg.generate_proof(coeffs, index_x)		#! Generating proof
	#print("proof =", proof)
	#tx = contract.verify(format_FQ_G1Point(commit),	#! Verifying on-chain
	#	    			format_FQ_G1Point(proof),
	#					format_field_to_int(index_x),
	#					format_field_to_int(value_y))
	#print("\nResult of on-chain verification:", tx)
	#print(kzg.verify(commit, proof, index_x, value_y))
def testNewton(x_values, values, index):
	kzg = KZG10()
	start = time.thread_time()
	coeffs = kzg.generate_coeffs_for2(values)
	#print("coeffs:", coeffs[0])
	y = kzg.newton_poly2(coeffs[0], x_values, index)
	end = time.process_time() - start
	print(f"Newton Interpolation for {len(values)} values")
	print('-'*50)
	print(f"Cpu time: {end} secs")
	print(f"Is evaluation correct? {int(y) == values[index]}")
	print('='*50)

def testSpline(x_values, values, index):
	kzg = KZG10()
	start = time.process_time()
	coeffs = kzg.cubic_spline_coefficients(x_values, values)
	y = kzg.evaluate_cubic_spline(coeffs, x_values, index)
	end = time.process_time() - start
	print(f"Cubic Spline Interpolation for {len(values)} values")
	print('-'*50)
	print(f"Cpu time: {end} secs")
	print(f"Is evaluation correct? {int(y) == values[index]}")

if __name__ == "__main__":
	#Prove()
	#testRun()
	#testMultiProof()
	values = [randint(1,500) for i in range(500)]
	x_values = [i for i in range(len(values))]
	testNewton(x_values, values, randint(0, len(values)-1))
	testSpline(x_values, values, randint(0, len(values)-1))



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

def get_initial_point(result: Point2D[Field], n: int) -> Point2D[Field]:
    while n > 1:
        if n % 2:
            result = add(result, curve.double(result))
        n = n // 2
        result = curve.double(result)
    return result

