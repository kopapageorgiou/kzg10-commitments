from typing import NamedTuple, Union, List
from functools import reduce
import operator

from KZG10 import PointG1, curve, Field, GF


class DivisionResult(NamedTuple):
    q: 'Polynomial'
    r: 'Polynomial'


class Polynomial(object):
    def __init__(self, F, coeffs=None):
        self.F = F
        coeffs = coeffs or []
        self.coeffs = [F(_) for _ in coeffs]

    def eval_g1(self, point_powers: List[PointG1]) -> PointG1:
        """
        Evaluate the polynomial at an abstract point
        Provide powers of the point, e.g.
            [g^0, g^1, g^2...]
        returns
            sum([(g^0) * c_0, g*c_1, (g^i)*c_i, ...])
        """
        assert len(point_powers) >= len(self.coeffs)
        result = None
        for i, c_i in enumerate(self.coeffs):
            x = point_powers[i]         
            y = curve.multiply(x, c_i)              
            result = y if result is None else curve.add(result, y)
        return result

    def eval_scalar(self, point: Field) -> Field:
        """
        Evaluate the polynomial at a known scalar point
        """
        return reduce(operator.add, [(point**i) * c_i for i, c_i in enumerate(self.coeffs)])

    def __repr__(self):
        x = [f'{int(c_i)}*x^{i}' for i, c_i in enumerate(self.coeffs) if i != 0]
        terms = ', '.join(x)
        return f'<Polynomial[{terms}]>'

    @classmethod
    def _compact(cls, coeffs) -> List[Field]:
        """Removes all trailing zero coefficients"""
        new_coeffs = list()
        skipping = True
        # TODO: use slicing instead of copying
        for c_i in coeffs[::-1]:
            if skipping:
                skipping = c_i == 0
            if not skipping:
                new_coeffs.append(c_i)
        return new_coeffs[::-1] 

    def compact(self) -> 'Polynomial':
        return Polynomial(self.F, self._compact(self.coeffs))

    def __call__(self, point: Union[List[PointG1], Field]) -> Union[PointG1, Field]:
        """
        Evaluate the polynomial at a specific point
        """
        if isinstance(point, list):
            assert isinstance(point[0], PointG1)
            return self.eval_g1(point)
        elif isinstance(point, Field):
            return self.eval_scalar(point)
        raise TypeError(f"Unknown how to evaluate at {type(point)}: {point}")

    def mul_poly(self, other: 'Polynomial') -> 'Polynomial':
        """
        From: https://jeremykun.com/tag/division-algorithm/
        """
        newCoeffs = [self.F(0) for _ in range(len(self) + len(other) - 1)] 
        for i,a in enumerate(self.coeffs):
            for j,b in enumerate(other.coeffs):
                newCoeffs[i+j] += a*b
        return Polynomial(self.F, newCoeffs).compact()

    def mul_scalar(self, other: Field) -> 'Polynomial':
        return Polynomial(self.F, [_ * other for _ in self.coeffs])

    def __mul__(self, other: Union[Field, 'Polynomial']) -> 'Polynomial':
        if isinstance(other, Polynomial):
            return self.mul_poly(other)
        elif isinstance(other, Field):
            return self.mul_scalar(other)
        raise TypeError(f"Unknown how to evaluate at {type(point)}: {point}")

    def __add__(self, other: 'Polynomial') -> 'Polynomial':
        new_length = max(len(self.coeffs), len(other.coeffs))
        result = [self.F(0)] * new_length
        for i in range(new_length):
            if len(self) > i:
                result[i] += self.coeffs[i]
            if len(other) > i:
                result[i] += other.coeffs[i]
        return Polynomial(self.F, result).compact()

    def __neg__(self) -> 'Polynomial':
        return Polynomial(self.F, [-c_i for c_i in self.coeffs])

    def __sub__(self, other: 'Polynomial') -> 'Polynomial':
        return self + -other

    def is_zero(self):
        return len(self.coeffs) == 0 or all([c_i == 0 for c_i in self.coeffs])

    @classmethod
    def random(cls, F: GF, n: int) -> 'Polynomial':
        return cls(F, [F.random() for _ in range(n)])

    def inverse(self):
        return Polynomial(self.F, [_.inverse() for _ in self.coeffs])

    def __iter__(self):
        return iter(self.coeffs)

    def __len__(self):
        return len(self.coeffs)

    def __getitem__(self, index: int) -> Field:
        return self.coeffs[index]

    def degree(self):
        return len(self) - 1

    def __eq__(self, other: 'Polynomial') -> 'Polynomial':
        return self.degree() == other.degree() and all([x==y for x, y in zip(self.coeffs, other.coeffs)])

    def leadingCoefficient(self) -> Field:
        return self.coeffs[-1]

    @classmethod
    def zero(cls, F):
        return cls(F, [])

    def __divmod__(self, other: 'Polynomial') -> 'Polynomial':
        # Synthetic division of polynomials
        # https://rosettacode.org/wiki/Polynomial_synthetic_division#Python
        """
        dividend = list(self.compact().coeffs[::-1])
        divisor = list(other.compact().coeffs[::-1])
        out = list(dividend)
        normalizer = divisor[0].inverse()
        for i in range(len(dividend) - (len(divisor)-1)):
            out[i] *= normalizer
            coef = out[i]
            if coef != 0:
                for j in range(1, len(divisor)):
                    out[i + j] += -divisor[j] * coef
        separator = -(len(divisor)-1)
        q = Polynomial(self.F, out[:separator][::-1]).compact()
        r = Polynomial(self.F, out[separator:][::-1]).compact()
        return DivisionResult(q, r)
        """
        # From: https://jeremykun.com/tag/division-algorithm/
        quotient, remainder = self.zero(self.F), self
        divisorDeg = other.degree()
        divisorLC = other.leadingCoefficient().inverse()

        while remainder.degree() >= divisorDeg:
            monomialExponent = remainder.degree() - divisorDeg
            monomialZeros = [self.F(0) for _ in range(monomialExponent)]
            monomialDivisor = Polynomial(self.F, monomialZeros + [remainder.leadingCoefficient() * divisorLC])

            quotient += monomialDivisor
            remainder -= monomialDivisor * other

        return DivisionResult(quotient, remainder)

    def __truediv__(self, other: 'Polynomial') -> 'Polynomial':
        q, r = divmod(self, other)
        return q


def test_division():
    Fp = GF(curve.curve_order)
    a = Polynomial(Fp, [14, 9, 1])
    b = Polynomial(Fp, [7, 1])
    q, r = divmod(a, b)
    assert q[0] == 2
    assert q[1] == 1
    assert r.is_zero() == True

    # (2x^2 - 5x - 1) / x-3
    a = Polynomial(Fp, [-1, -5, 2])
    b = Polynomial(Fp, [-3, 1])
    q, r = divmod(a, b)
    assert len(q) == 2
    assert len(r) == 1
    assert q[0] == 1
    assert q[1] == 2
    assert r[0] == 2

    # sage: f(x) = 18*x^3 - x^2 - 5*x - 8
    # sage: g(x) = x^2 - 3
    # sage: f.maxima_methods().divide(g)
    # [18*x - 1, 49*x - 11] ... ?
    # Instead, just check quotient remainder theorem
    a = Polynomial(Fp, [-8, -5, -1, 18])
    b = Polynomial(Fp, [3, 1])
    q, r = divmod(a, b)
    c = (b*q) + r
    assert a == c


def test_multiplication():
    Fp = GF(curve.curve_order)
    a = Polynomial(Fp, [Fp(_) for _ in [-5, 1]])
    b = Polynomial(Fp, [Fp(_) for _ in [2, 1]])
    c = a * b
    assert len(c) == 3
    assert c[0] == Fp(-10)
    assert c[1] == Fp(-3)
    assert c[2] == Fp(1)


def main():
    test_division()
    test_multiplication()


if __name__ == "__main__":
    main()