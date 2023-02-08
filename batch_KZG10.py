from functools import reduce
import operator
from BlockchainLib import *
from BLS import *
from KZG10 import TrustedSetup, CommitDivision, curve, GF, polynomial, CommitSum


"""
The batch commitment scheme,
using a `gamma` variable randomly chosen by the verifier
We then multiply the commitments by the powers of gamma to 
commit to multiple polynomials.
Say, for example we have two polynomials with 3 coefficients each
sage:
phi = lambda k, poly: poly[0] + ((k^1)*poly[1]) + ((k^2)*poly[2])
psi = lambda a, b, poly: ((phi(a, poly) - phi(b, poly)) / (a - b))
poly0 = var(['p0c%d' % _ for _ in range(3)])
poly1 = var(['p1c%d' % _ for _ in range(3)])
poly2 = var(['p2c%d' % _ for _ in range(3)])
x, z, gamma = var('x z gamma')
# Then make witnesses for the two polynomials evaluated at the same point 
h_0 = (gamma^1) * psi(x, z, poly0)
h_1 = (gamma^2) * psi(x, z, poly1)
h_2 = (gamma^3) * psi(x, z, poly2)
W = H = h_0 + h_1 + h_2
# Where `H` is equivalent to the batch witness `W`
# We can see the form the batch witness polynomial takes
sage: factor(h_0)
(p0c2*x + p0c2*z + p0c1)*gamma
sage: factor(h_0 + h_1)
(gamma*p1c2*x + gamma*p1c2*z + gamma*p1c1 + p0c2*x + p0c2*z + p0c1)*gamma
sage: factor(h_0 + h_1 + h_2)
(gamma^2*p2c2*x + gamma^2*p2c2*z + gamma^2*p2c1 + gamma*p1c2*x + gamma*p1c2*z + gamma*p1c1 + p0c2*x + p0c2*z + p0c1)*gamma
# Commit to the polynomials
cm_0 = phi(x, poly0)
cm_1 = phi(x, poly1)
cm_2 = phi(x, poly2)
# Then verifier computes the elements F and v
F_0 = (gamma^1) * cm_0
F_1 = (gamma^2) * cm_1
F_2 = (gamma^3) * cm_2
F = F_0 + F_1 + F_2
# For each of the commitments
(gamma^1) * (cm_0 - phi(z, poly0)) == h_0 * (x-z)
(gamma^2) * (cm_1 - phi(z, poly1)) == h_1 * (x-z)
(gamma^3) * (cm_2 - phi(z, poly2)) == h_2 * (x-z)
# Symbolic computations of the evaluations... 
s0 = gamma^1 * phi(z, poly0)
s1 = gamma^2 * phi(z, poly1)
s2 = gamma^3 * phi(z, poly2)
v = s0 + s1 + s2
# The pairing equation verifies that
# F-v == W * (x-z)
# Looking at the expansion, we can see they are equivalent
sage: factor(expand(W * (x-z)))
(gamma^2*p2c2*x + gamma^2*p2c2*z + gamma^2*p2c1 + gamma*p1c2*x + gamma*p1c2*z + gamma*p1c1 + p0c2*x + p0c2*z + p0c1)*gamma*(x - z)
sage: factor(expand(F-v))
(gamma^2*p2c2*x + gamma^2*p2c2*z + gamma^2*p2c1 + gamma*p1c2*x + gamma*p1c2*z + gamma*p1c1 + p0c2*x + p0c2*z + p0c1)*gamma*(x - z)
sage: factor((F-v) / (W * (x-z)))
1
"""

def main():
	F = GF(curve.curve_order)
	n_polynomials = 3
	n_coeffs = 3

	poly_coeffs = [ [F.random() for _ in range(n_coeffs)] for i in range(n_polynomials) ]
	print(poly_coeffs)
	PK = TrustedSetup.generate(F, n_coeffs, True)

	# Commit to polynomials, phi(f_i, X)
	cm = []
	for coeffs_i in poly_coeffs:
		cm.append( CommitSum(PK, coeffs_i) )

	# Prover creates a random scalar
	z = F.random()
	print("z", z)
	# Verifier sends random scalar
	gamma = F.random()
	print("gamma", gamma)
	# Prover computes the polynomial commitment h(X)
	s = []
	W = None
	for i, coeffs_i in enumerate(poly_coeffs):
		s.append(polynomial(z, coeffs_i))  # s_i = f_i(z)
		h_i = CommitDivision(PK, z, coeffs_i) # h_i = psi_z(f_i, X)
		g1_h_i = curve.multiply(h_i, int(gamma**(i+1))) # y^i * h_i
		W = g1_h_i if W is None else curve.add(W, g1_h_i)
	# Prover sends W, z and s to verifier
	print("W", W)
	# Verifier computes the element F
	g1_F = None
	for i, cm_i in enumerate(cm):
		F_i = curve.multiply(cm_i, int(gamma**(i+1)))
		g1_F = F_i if g1_F is None else curve.add(g1_F, F_i)
	print("g1_F", g1_F)
	# And verifier computes the element `v`
	v = reduce(operator.add, [(gamma**(i+1)) * s_i for i, s_i in enumerate(s)])
	print("v", v)
	g1_v = curve.multiply(curve.G1, int(v))
	print("g1_v", g1_v)

	# Then performs the pairing equation
	g1_F_minus_v = curve.add(g1_F, curve.neg(g1_v))
	print("g1_F_minus_v", g1_F_minus_v)
	g2_z = curve.multiply(curve.G2, int(z))
	g2_x_minus_z = curve.add(PK.g2_powers[1], curve.neg(g2_z))
	a = curve.pairing(curve.G2, g1_F_minus_v)
	b = curve.pairing(g2_x_minus_z, W)
	print('a', a)
	print('b', b)
	c = curve.pairing(curve.neg(g2_x_minus_z), W)
	print('a*c', a*c)


if __name__ == "__main__":
	main()