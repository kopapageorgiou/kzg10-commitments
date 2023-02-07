from typing import List, NamedTuple, Tuple
from functools import reduce
import operator

from KZG10 import TrustedSetup, CommitDivision, curve, GF, polynomial, CommitSum, Field, PointG1


"""
As per the batch version of KZG10
This opens many polynomials at many points
sage:
phi = lambda k, poly: poly[0] + ((k^1)*poly[1]) + ((k^2)*poly[2])
psi = lambda a, b, poly: ((phi(a, poly) - phi(b, poly)) / (a - b))
# Proving two polynomials, opened at different positions
polynomials = [var(['p%dc%d' % (_, j) for j in range(3)])
			   for _ in range(2)]
x = var('x')
z1, z2 = var('z1 z2')
s1 = phi(z1, polynomials[0])
s2 = phi(z2, polynomials[1])
cm1 = phi(x, polynomials[0])
cm2 = phi(x, polynomials[1])
gamma1, gamma2 = var('gamma1 gamma2')
W1 = gamma1 * psi(x, z1, polynomials[0])
W2 = gamma2 * psi(x, z2, polynomials[1])
r1, r2 = var('r1 r2')
v1 = (gamma1*s1)
v2 = (gamma2*s2)
F_1 = (gamma1 * cm1) - v1
F_2 = (gamma2 * cm2) - v2
F = (r1*F_1) + (r2*F_2)
lhs = F + ((r1*z1)*W1) + ((r2*z2)*W2)
rhs = ((r1 * W1) + (r2*W2)) * x
sage: factor((r1*W1+r2*W2)*x)
(gamma1*p0c2*r1*x + gamma2*p1c2*r2*x + gamma1*p0c2*r1*z1 + gamma2*p1c2*r2*z2 + gamma1*p0c1*r1 + gamma2*p1c1*r2)*x
sage: factor(expand(lhs))
(gamma1*p0c2*r1*x + gamma2*p1c2*r2*x + gamma1*p0c2*r1*z1 + gamma2*p1c2*r2*z2 + gamma1*p0c1*r1 + gamma2*p1c1*r2)*x
"""

Coefficients = List[Field]
PolynomialGroup = List[Coefficients]


class Commitment(NamedTuple):
	c: PolynomialGroup
	cm: List[PointG1]
	z: Field
	s: List[Field]

	@classmethod
	def create(cls, PK: TrustedSetup, polysys: PolynomialGroup):
		z = PK.F.random()
		s = [polynomial(z, _) for _ in polysys]
		cm = [CommitSum(PK, _) for _ in polysys]
		return cls(polysys, cm, z, s)

	def witness(self, PK: TrustedSetup, gamma: Field) -> PointG1:
		return reduce(curve.add, [curve.multiply(CommitDivision(PK, self.z, coeffs_i), int(gamma**(i+1)))
								  for i, coeffs_i in enumerate(self.c)])

	def prepare_F_and_v(self, gamma: Field) -> Tuple[PointG1, Field]:
		g1_F = reduce(curve.add, [curve.multiply(cm_i, int(gamma**(i+1))) for i, cm_i in enumerate(self.cm)])
		v = reduce(operator.add, [(gamma**(j+1)) * s_i for j, s_i in enumerate(self.s)])
		return g1_F, v

	def prepare_F_minus_v(self, gamma: Field) -> PointG1:
		g1_F, v = self.prepare_F_and_v(gamma)
		g1_v = curve.multiply(curve.G1, int(v))
		return curve.add(g1_F, curve.neg(g1_v))


def Step1_Prover(PK: TrustedSetup, many_polynomials: List[PolynomialGroup]) -> List[Commitment]:
	"""
	Commit to many groups of polynomials
	"""
	return [Commitment.create(PK, _) for _ in many_polynomials]


def Step2_Verifier(PK: TrustedSetup, step1: List[Commitment]) -> List[Field]:
	"""
	Verifier challenges prover, by providing random gamma values for each commitment
	"""
	return [PK.F.random() for _ in range(len(step1))]


def Step3_Prover(PK: TrustedSetup, step1: List[Commitment], step2: List[Field]) -> List[PointG1]:
	"""
	Prover commits to the groups of polynomials, using the gamma values from the verifier in step2
	"""
	assert len(step1) == len(step2)
	return [_.witness(PK, gamma) for gamma, _ in zip(step2, step1)]


def Step4_Verifier_Individually(PK: TrustedSetup, step1: List[Commitment], step2: List[Field], step3: List[PointG1]):
	results = []
	for group, gamma, W in zip(step1, step2, step3):
		g1_F_minus_v = group.prepare_F_minus_v(gamma)
		g2_z = curve.multiply(curve.G2, int(group.z))
		g2_x_minus_z = curve.add(PK.g2_powers[1], curve.neg(g2_z))
		a = curve.pairing(curve.G2, g1_F_minus_v)
		b = curve.pairing(curve.neg(g2_x_minus_z), W)
		results.append(a*b == curve.FQ12.one())
	return all(results)


def Step4_Verifier(PK: TrustedSetup, step1: List[Commitment], step2: List[Field], step3: List[PointG1]):
	lhs_args = list()
	rhs_args = list()
	effs = list()
	for group, gamma, W in zip(step1, step2, step3):
		r = PK.F.random()
		g1_F_minus_v = group.prepare_F_minus_v(gamma)
		effs.append(curve.multiply(g1_F_minus_v, int(r)))
		lhs_args.append(curve.multiply(W, int(r*group.z)))
		rhs_args.append(curve.multiply(W, int(r)))

	lhs_args.insert(0, reduce(curve.add, effs))

	lhs = reduce(curve.add, lhs_args)
	rhs = reduce(curve.add, rhs_args)

	a = curve.pairing(curve.G2, lhs)
	b = curve.pairing(PK.g2_powers[1], curve.neg(rhs))
	print('a', a)
	print('b', b)
	return a*b == curve.FQ12.one()


def main():
	F = GF(curve.curve_order)
	n_polynomials = 3
	n_coeffs = 3
	n_openings = 2

	groups = [[[F.random() for _ in range(n_coeffs)]
			    for i in range(n_polynomials)]
			    for j in range(n_openings)]

	PK = TrustedSetup.generate(F, n_coeffs, True)
	step1 = Step1_Prover(PK, groups)
	step2 = Step2_Verifier(PK, step1)
	step3 = Step3_Prover(PK, step1, step2)

	print('Verifying each opening individually')
	print(Step4_Verifier_Individually(PK, step1, step2, step3))

	print('Verifying all openings together in a batch')
	step4 = Step4_Verifier(PK, step1, step2, step3)
	print(step4)


if __name__ == "__main__":
	main()