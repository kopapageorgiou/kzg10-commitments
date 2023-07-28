from typing import List, Tuple, Type
import json
from functools import reduce
from py_ecc.fields import bn128_FQ as FQ
from py_ecc.fields import bn128_FQ2 as FQ2
from py_ecc.fields import bn128_FQ12 as FQ12
from py_ecc.bn128 import bn128_pairing as curvePairing
from py_ecc import bn128 as curve
from py_ecc.bn128 import bn128_curve as cu
from py_ecc.bn128.bn128_curve import G2, G1
from py_ecc.typing import Field
from random import randint
import numpy as np

G1_POINTS = []
G2_POINTS = []

SRS_G2_1 = (FQ2([
                0x04c5e74c85a87f008a2feb4b5c8a1e7f9ba9d8eb40eb02e70139c89fb1c505a9,
                0x21a808dad5c50720fb7294745cf4c87812ce0ea76baa7df4e922615d1388f25a
            ]),
            FQ2([
                0x2d58022915fc6bc90e036e858fbc98055084ac7aff98ccceb0e3fde64bc1a084,
                0x204b66d8e1fadc307c35187a6b813be0b46ba1cd720cd1c4ee5f68d13036b4ba
                
            ]))
class KZG10(object):

    NEWTON = "newton"
    LAGRANGE = "lagrange"
    MONOMIAL = "monomial"

    def __init__(self, field = cu.curve_order) -> None:
        """
        Load G1 Points from tau json file
        and create GF with field= bn128.curve_order a.k.a BABYJUB_P (the same field with the one in solidity)
        """
        with open('taug1_65536.json', 'r') as file:
            g1_points = json.load(file)
        file.close()
        
        for point in g1_points:
            G1_POINTS.append(format_to_FQ(int(point[0], base=16), int(point[1], base=16)))
        
        with open('taug2_65536.json', 'r') as file:
            g2_points = json.load(file)
        file.close()

        for i in range(len(g2_points)):
            G2_POINTS.append((FQ2([int(g2_points[i][0], base=16), int(g2_points[i][1], base=16)]),
                             FQ2([int(g2_points[i][2], base=16), int(g2_points[i][3], base=16)])))
        
        self.field = field
        self.F = GF(field)

    def choose_method(self, method: str):
        if method == self.NEWTON:
            return Newton(self.field)
        elif method == self.LAGRANGE:
            return LaGrange(self.field)
        elif method == self.MONOMIAL:
            return Monomial(self.field)
        
    def evalPolyAt(self, coefficients: List[Field], index: Field):
        result = self.F(0)
        power_of_x = self.F(1)

        for coeff in coefficients:
            result = (result + (power_of_x * coeff))
            power_of_x = (power_of_x * index)

        return result
        """result = self.F(0)
        n = len(coefficients)
        for i in range(n):
            result += coefficients[i] * (index ** (n-i-1))"""

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
        return _swap(pol)
    
    def newton_interpolation(self, values: List[int]):
        x_values = []
        val = []
        for x, y in enumerate(values):
            x_values.append(x)
            val.append(y)
        
        pol = self.divided_diff([self.F(x) for x in x_values], [self.F(y) for y in val])[0]
        return pol

    def divided_diff(self, x,y):
        n = len(y)

        # coef = np.zeros([n,n])
        #coef = [[self.F(0)]]
        coef = [[self.F(0) for i in range(n)] for j in range(n)]
        
        # coef[:,0] = y
        # first_col = []
        # for row in coef:
        #     first_col.append(row[0])
        for i in range(n):
            coef[i][0] = Field(y[i],self.field)

        for j in range(1,n):
            for i in range(n-j):
                coef[i][j] = (coef[i+1][j-1]-coef[i][j-1])/((x[i+j]-x[i]))
        return coef
        
    def newton_poly(self, coef, x):
        x_data = [i for i in range(len(coef))]
        n = len(x_data) - 1
        p = coef[n]

        for k in range(1, n+1):
            p = coef[n-k] + (x-x_data[n-k])*p
        return p
    
    def eval_newton_poly_at(self, coef, x):
        x_data = [i for i in range(len(coef))]
        n = len(x_data) - 1
        p = coef[n]

        for k in range(1, n+1):
            p = coef[n-k] + (self.F(x)-self.F(x_data[n-k]))*p
            
        return p
    

    def generate_proof(self, coefficients: List[Field], index: Field):
        quotientCoefficients = self._genQuotientPolynomial(coefficients, index)
        
        return self.generate_commitment(quotientCoefficients)
        #return self.custom_commit(quotientCoefficients)

    def _genQuotientPolynomial(self, coefficients: List[Field], xVal: Field):
        yVal = self.evalPolyAt(coefficients, xVal)
        print("Y-val quot:", yVal)
        x = [self.F(0), self.F(1)]
        res = self._divPoly(self._subPoly(coefficients, [yVal]), self._subPoly(x, [xVal]))[0]
        return res
    
    def get_index_x(self, index: int):
        return self.F(index)
    
    def custom_generate_proof(self, coefficients: List[Field], index: Field):
        quotientCoefficients = self._customgenQuotientPolynomial(coefficients, index)
        #print("here", quotientCoefficients)
        return self.generate_commitment(quotientCoefficients)

    def _customgenQuotientPolynomial(self, coefficients: List[Field], xVal: int):
        yVal = self.evaluate_spline_point(coefficients, xVal)
        
        #segment_coefficients = self.pre_coeffs(coefficients, xVal)

        #yVal = self.evaluate_cubic_spline_seg(segment_coefficients, xVal)
        print("Y-val quot:", yVal)
        x = [self.F(0), self.F(1)]
        #print("Segment_coefficients:", segment_coefficients)
        res = self._divPoly(self._subPoly(coefficients, [yVal]), self._subPoly(x, [xVal]))[0] # type: ignore
        print("Res", res)
        return res
    
    def pre_coeffs(self, coeffs, xVal):
        n = len(coeffs)
        x_data = [i for i in range(n)]
        for i in range(n-1):
            if x_data[i] <= xVal <= x_data[i+1]:
                break
        
        coefficients = []
        coefficients.append(coeffs[4*i])
        coefficients.append(coeffs[4*i+1])
        coefficients.append(coeffs[4*i+2])
        coefficients.append(coeffs[4*i+3])
        #print(type(h), type(t))
        return coefficients


    """
	Copute commitment to the evaluation of the polynomial given the coefficients
	"""
    def generate_commitment(self, coefficients: List[Field]):
        return reduce(cu.add, [cu.multiply(G1_POINTS[i], int(c_i))
							        for i, c_i in enumerate(coefficients)])
    def custom_commit(self, coeffs):
        return _reduce(_add, [_multiply(G1_POINTS[i], int(c_i))
                        for i, c_i in enumerate(coeffs)])
    
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
        # Determine the degrees of the input polynomials
        degree1 = len(poly1) - 1
        degree2 = len(poly2) - 1

        # Determine the degree of the resulting polynomial
        result_degree = degree1 + degree2

        # Initialize the result polynomial with zeros
        result_poly = [self.F(0)] * (result_degree + 1)

        # Perform polynomial multiplication
        for i in range(degree1 + 1):
            for j in range(degree2 + 1):
                result_poly[i + j] += poly1[i] * poly2[j]

        return result_poly
    
    def div_polys(self, a, b):
        """
        Long polynomial difivion for two polynomials in coefficient form
        """
        a = [x for x in a]
        o = []
        apos = len(a) - 1
        bpos = len(b) - 1
        diff = apos - bpos
        while diff >= 0:
            quot = a[apos] / self.F(b[bpos])
            o.insert(0, quot)
            for i in range(bpos, -1, -1):
                a[diff + i] -= self.F(b[i]) * quot
            apos -= 1
            diff -= 1
        return o

    def _mulPoly(self, poly1, poly2):
    # Initialize an output polynomial with all coefficients set to zero
        output_poly = [self.F(0)] * (len(poly1) + len(poly2) - 1)
        # Multiply each term of poly1 with each term of poly2, and add the results to the output polynomial
        for i, coeff1 in enumerate(poly1):
            for j, coeff2 in enumerate(poly2):
                #print(coeff1, coeff2)
                output_poly[i+j] += coeff1 * coeff2

        return output_poly


    """def _lagrange_interpolation(self, x_vals, y_vals):
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
        return coefficients"""
    
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
    

    def generate_multi_proof(self, coefficients: List[Field], indices: List[int], values: List[int]):
        indices = [self.F(indice) for indice in indices]
        values = [self.F(value) for value in values]
        poly = coefficients
        ipoly = self._genInterpolatingPoly(poly, indices)
        #print(ipoly)
        zpoly = self._genZeroPoly(indices)
        #print(zpoly)
        qPoly = self._divPoly(self._subPoly(poly, ipoly), zpoly)
        #print(qPoly)
        multiproof = self._polyCommit(qPoly[0])
        iCoeffs = _swap(self._lagrange_inter(indices, values))
        #print("here",iCoeffs)
        #print(ipoly)
        return multiproof, iCoeffs, zpoly

    
    def _polyCommit(self, coefficients: List[Field]):
        return reduce(cu.add, [cu.multiply(G2_POINTS[i], int(c_i))
							        for i, c_i in enumerate(coefficients)])

    def _srsG2(self, amount: int):
        return [G2_POINTS[i] for i in range(amount)]
    
    def _genZeroPoly(self, indices: List[Field]):
        zPoly = [self.F(-1) * indices[0], self.F(1)]
        #print(zPoly)
        for indice in indices:
            zPoly = self._mulPoly(zPoly, [self.F(-1) * indice, self.F(1)])
            #print(zPoly[1:])

        return zPoly[1:]



    def mod(self, value: int) -> int:
        if value >= 0:
            return self.F(value % self.field)
        else:
            return self.F(((value % self.field) + self.field) % self.field)

    def _genInterpolatingPoly(self, poly: List[Field], indices: List[Field]):
        x = []
        values = []
        
        for indice in indices:
            values.append(self.evalPolyAt(poly, indice))
            x.append(indice)
        
        coeffs = _swap(self._lagrange_inter(x, values))
        #print(self.evalPolyAt(coeffs, self.F(4)))
        return coeffs
            
    def verify_off_chain(self, commitment, proof, index, value):
        commitmentMinusA = cu.add(commitment, cu.neg(cu.multiply(cu.G1, int(value))))
        negProof = cu.neg(proof)
        indexMulProof = cu.multiply(proof, int(index))
        #print("here",indexMulProof)
        #return [commitmentMinusA, negProof, indexMulProof]
        a = curvePairing.pairing(G2, cu.add(indexMulProof, commitmentMinusA))
        b = curvePairing.pairing(SRS_G2_1, negProof)
        ab = a*b
        #print(a)
        #print(b)
        #print("ab", ab)
        return ab == FQ12.one()
    
    def cubic_spline_coefficients(self, x, y):
        n = len(x)
        h = [self.F(0)] * (n-1)
        for i in range(n-1):
            h[i] = self.F(x[i+1]) - self.F(x[i])

        alpha = [self.F(0)] * (n-1)
        for i in range(1, n-1):
            alpha[i] = (self.F(3)/h[i]) * (self.F(y[i+1])-self.F(y[i])) - (self.F(3)/h[i-1]) * (self.F(y[i])-self.F(y[i-1]))

        l = [self.F(1)] + [self.F(0)] * (n-2)
        mu = [self.F(0)] * (n-1)
        z = [self.F(0)] * n
        for i in range(1, n-1):
            l[i] = self.F(2) * (self.F(x[i+1]) - self.F(x[i-1])) - h[i-1] * mu[i-1]
            mu[i] = h[i] / l[i]
            z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i]

        c = [self.F(0)] * n
        b = [self.F(0)] * (n-1)
        d = [self.F(0)] * (n-1)
        c[-1] = self.F(0)
        for j in range(n-2, -1, -1):
            c[j] = z[j] - mu[j]*c[j+1]
            b[j] = (self.F(y[j+1]) - self.F(y[j])) / h[j] - h[j] * (c[j+1] + self.F(2)*c[j]) / self.F(3)
            d[j] = (c[j+1] - c[j]) / (self.F(3)*h[j])

        # The coefficients of the cubic polynomials are stored in a list
        coeffs = []
        for i in range(n-1):
            coeffs.append(self.F(y[i]))
            coeffs.append(b[i])
            coeffs.append(c[i])
            coeffs.append(d[i])
        return coeffs
    
    def generate_spline_monomial_coeffs(self, x, y):
        coeffs = self.cubic_spline_coefficients(x, y)
        num_segments = len(x) - 1

        monomial_coeffs = []
        for i in range(num_segments):
            y_i, b_i, c_i, d_i = coeffs[i * 4 : (i + 1) * 4]

            # Generate the monomial coefficients for the cubic polynomial on this segment
            monomial_coeffs.append(y_i)
            monomial_coeffs.append(b_i)
            monomial_coeffs.append(c_i / self.F(2))
            monomial_coeffs.append(d_i / self.F(3))

        return monomial_coeffs

    def evaluate_spline_point(self, monomial_coeffs, xVal):
        n = len(monomial_coeffs) // 4  # Number of segments in the spline
        x = [i for i in range(n)]
        num_segments = len(x) - 1

        # Find the segment index that contains xVal
        segment_index = 0
        for i in range(num_segments):
            if x[i] <= xVal <= x[i + 1]:
                segment_index = i
                break
        #xVal = self.F(xVal)
        #x = [self.F(xv) for xv in x]
        # Extract the monomial coefficients for the segment
        a, b, c, d = monomial_coeffs[segment_index * 4 : (segment_index + 1) * 4]
        x_i = x[segment_index]

        # Evaluate the cubic polynomial at xVal
        dx = xVal - x_i
        result = a + b * dx + (c * self.F(2)) * dx**2 + (d * self.F(3)) * dx**3

        return result
    
    def generate_spline_polynomial(self, x, y):
        monomial_coeffs = self.generate_spline_monomial_coeffs(x, y)
        num_segments = len(x) - 1

        # Initialize the polynomial coefficients
        poly_coeffs = [self.F(0)] * (len(x)*4)

        # Combine the monomial coefficients for each segment to form the polynomial
        for i in range(num_segments):
            a, b, c, d = monomial_coeffs[i * 4 : (i + 1) * 4]

            # Accumulate the contributions from each segment
            poly_coeffs[i] += a
            poly_coeffs[i + 1] += b
            poly_coeffs[i + 2] += c / self.F(2)
            poly_coeffs[i + 3] += d / self.F(3)

        # The resulting poly_coeffs list now represents the polynomial
        return poly_coeffs[:len(x)]




    def evaluate_cubic_spline_seg(self, coeffs, x):
        """
        Evaluates the cubic spline polynomial at a given value of x, using the coefficients
        generated by the cubic_spline_coefficients function.

        :param coeffs: A list of coefficients for a specific segment generated by the cubic_spline_coefficients function.
        :param x: The value of x at which to evaluate the cubic spline.
        :return: The interpolated value of f(x) using the cubic spline polynomial.
        """
        x = x % 4
        if len(coeffs) != 4:
            raise ValueError("Invalid number of coefficients for cubic spline interpolation.")

        a = coeffs[0]
        b = coeffs[1]
        c = coeffs[2]
        d = coeffs[3]

        h = self.F(1)  # Assuming the x data for each segment is normalized to [0, 1]
        #if x==0:
        t = (self.F(x) - self.F(0)) / h
        #else:
        #    t = (self.F(x) - self.F(x-1)) / h

        return a + b * t + c * t ** 2 + d * t ** 3
    
    def evaluate_cubic_spline(self, coeffs, x):
        """
        Evaluates the cubic spline polynomial at a given value of x, using the coefficients
        generated by the cubic_spline_coefficients function.
        
        :param x: The value of x at which to evaluate the cubic spline.
        :param coeffs: A list of coefficients generated by the cubic_spline_coefficients function.
        :param x_data: A list of the x values used to generate the coefficients.
        :return: The interpolated value of f(x) using the cubic spline polynomial.
        """
        
        # Find the interval containing x
        n = len(coeffs)
        x_data = [i for i in range(n)]
        for i in range(n-1):
            if x_data[i] <= x <= x_data[i+1]:
                break
        
        # Evaluate the polynomial at x
        h = self.F(x_data[i+1]) - self.F(x_data[i])
        t = (self.F(x) - self.F(x_data[i])) / h
        
        a = coeffs[4*i]
        b = coeffs[4*i+1]
        c = coeffs[4*i+2]
        d = coeffs[4*i+3]
        #print(type(h), type(t))
        return a + b*t + c*t**2 + d*t**3
    
    def b_spline_interpolation(self, x_data, y_data, degree):
        n = len(x_data)
        m = n + degree + 1

        # Calculate the knot vector (using chord-length knots in this case)
        t = [self.F(0)] * m
        for i in range(1, m):
            t[i] = t[i - 1] + ((self.F(x_data[n - 1]) - self.F(x_data[0])) / (self.F(m) - self.F(1)))

        # Set up the linear system
        A = [[self.F(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                A[i][j] = self.basis_function(t, j, degree, x_data[i])

        # Solve the linear system
        coefficients = self.gauss_elimination(A, y_data)

        return coefficients

    def basis_function(self, t, i, k, x):
        if k == 0:
            return self.F(1) if int(t[i]) <= x < int(t[i + 1]) else self.F(0)

        if t[i + k] == t[i]:
            c1 = self.F(0)
        else:
            c1 = (self.F(x) - t[i]) / (t[i + k] - t[i]) * self.basis_function(t, i, k - 1, x)

        if t[i + k + 1] == t[i + 1]:
            c2 = self.F(0)
        else:
            c2 = (t[i + k + 1] - self.F(x)) / (t[i + k + 1] - t[i + 1]) * self.basis_function(t, i + 1, k - 1, x)

        return c1 + c2


    def gauss_elimination(self, A, b):
        n = len(A)
        b = [self.F(v) for v in b]
        # Forward elimination
        for i in range(n):
            pivot = A[i][i]
            for j in range(i + 1, n):
                factor = A[j][i] / pivot
                for k in range(i, n):
                    A[j][k] -= factor * A[i][k]
                b[j] -= factor * b[i]

        # Backward substitution
        x = [self.F(0)] * n
        for i in range(n - 1, -1, -1):
            x[i] = b[i]
            for j in range(i + 1, n):
                x[i] -= A[i][j] * x[j]
            x[i] /= A[i][i]

        return x


    def linear_inter_evaluation(self, coeffs, x_data, index):
        x_data = [self.F(x) for x in x_data]
        """
        Given the coefficients of a linear interpolation function, the x values of the input points,
        and an index of the point to evaluate, returns the y value of the interpolated point at that index.
        """
        if index < 0 or index > len(x_data)-1:
            raise ValueError("Index out of range.")
        
        slope, y_intercept = coeffs[index]
        x = x_data[index]
        return slope * x + y_intercept
   
    def linear_interpolation(self, x_values, y_values):
        x_values = [self.F(x) for x in x_values]
        y_values = [self.F(y) for y in y_values]
        """
        Given a list of points, returns a list of coefficients for a linear interpolation function.
        Each point in the input list should be a tuple of two numbers, representing x and y values respectively.
        """
        n = len(x_values)
        if n < 2:
            raise ValueError("At least two points are required.")

        # Calculate slopes and y-intercepts between each pair of adjacent points
        slopes = [(y_values[i+1] - y_values[i]) / (x_values[i+1] - x_values[i]) for i in range(n-1)]
        y_intercepts = [y_values[i] - slopes[i]*x_values[i] for i in range(n-1)]

        # Return a list of coefficients as tuples (slope, y-intercept)
        return [(slopes[i], y_intercepts[i]) for i in range(n-1)]



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


"""
* Custom integration of py_ecc and functools methods
"""
# Check if a point is the point at infinity
def _is_inf(pt):
    return pt is None

# Elliptic curve doubling
def _double(pt):
    if _is_inf(pt):
        return pt
    x, y = pt
    m = 3 * x**2 / (2 * y)
    newx = m**2 - 2 * x
    newy = -m * newx + m * x - y
    return (newx, newy)

# Elliptic curve point multiplication
def _multiply(pt, n: int):
    if n == 0:
        return None
    elif n == 1:
        return pt
    elif not n % 2:
        return _multiply(_double(pt), n // 2)
    else:
        return _add(_multiply(_double(pt), int(n // 2)), pt)
    
# Elliptic curve addition
def _add(p1, p2):
    if p1 is None or p2 is None:
        return p1 if p2 is None else p2
    x1, y1 = p1
    x2, y2 = p2
    if x2 == x1 and y2 == y1:
        return _double(p1)
    elif x2 == x1:
        return None
    else:
        m = (y2 - y1) / (x2 - x1)
    newx = m**2 - x1 - x2
    newy = -m * newx + m * x1 - y1
    assert newy == (-m * newx + m * x2 - y2)
    return (newx, newy)

def _reduce(func, data: List):
    x = data[0]

    for i in range(1,len(data)):
        x = func(x,data[i])
    return x

def format_to_FQ(x_point: int, y_point: int):
    return (FQ(x_point), FQ(y_point))

def format_FQ_G1Point(data: Tuple[Field, Field]):
    x, y = data
    return (int(str(x)), int(str(y)))

def format_proof(data):
    p1, p2 = data
    x1, x2 = p1.coeffs
    y1, y2 = p2.coeffs
    return ([int(x2),int(x1)],[int(y2),int(y1)])

def format_field_to_int(value: Field | List[Field]):
    if (type(value) == Field):
        return int(value) # type: ignore
    else:
        return [int(val) for val in value] # type: ignore
    
def _swap(l):
    def swapPositions(list, pos):
        #print(list[pos], list[pos-1])
        list[pos], list[-pos-1] = list[-pos-1], list[pos]
        return list
    res = []
    if len(l) % 2 == 1:
        middle = len(l)+1//2
        #print(middle)
        for i in range(middle-2):
            res = swapPositions(l, i).copy()
    else:
        middle = len(l)//2
        for i in range(middle):
            res = swapPositions(l,i).copy()
    return res
    

def main():
    kzg = KZG10(127)
    coeffs = kzg.generate_coeffs_for([5, 25, 125])
    commit = kzg.generate_commitment(coeffs)
    print(coeffs)
    for i in range(len(coeffs)):
        print(kzg.evalPolyAt(coeffs, kzg.F(i)))
    x = kzg.get_index_x(1)
    y = kzg.evalPolyAt(coeffs, x)

    proof = kzg.generate_proof(coeffs, x)
    kzg.verify_off_chain(commit, proof, x, y)

#** Currently Working
class Newton(KZG10):
    def __init__(self, field=cu.curve_order) -> None:
        self.field = field
        self.F = GF(field)
    
    def interpolate(self, values: List[int]):
        x_values = []
        val = []
        for x, y in enumerate(values):
            x_values.append(x)
            val.append(y)
        x_values = [self.F(x) for x in x_values]
        pol = super().divided_diff([self.F(x) for x in x_values], [self.F(y) for y in val])[0]
        newton = []
        for ck, xk in zip(pol[::-1], x_values[::-1]):
            newton = super()._addPoly(super()._mulPoly(newton, [-xk, self.F(1)]), [ck])
        #print(newton)
        return newton
    
    def eval_poly_at(self, coeffs: List[Field], index: Field):
        result = self.F(0)
        for i, coeff in enumerate(coeffs):
            result += coeff * (index ** i)
        return result

    def generate_commitment(self, coefficients):
        return super().generate_commitment(coefficients)
    
    def generate_proof(self, coefficients: List[Field], index: Field):
        quotientCoefficients = self._generate_quotient_polynomial(coefficients, index)
        #print("here", quotientCoefficients)
        return super().generate_commitment(quotientCoefficients)

    def _generate_quotient_polynomial(self, coefficients: List[Field], xVal: Field):
        yVal = self.eval_poly_at(coefficients, xVal)
        #print(xVal)
        #print("Y-val quot:", yVal)
        x = [self.F(0), self.F(1)]
        #print(super().div_polys(coefficients, [-xVal, 1]))
        res = super()._divPoly(super()._subPoly(coefficients, [yVal]), super()._subPoly(x, [xVal]))[0] # type: ignore
        #print("Res", res)
        return res


#** Currently working
class LaGrange(KZG10):
    def __init__(self, field=cu.curve_order) -> None:
        self.field = field
        self.F = GF(field)

    def interpolate(self, values: List[int]):
        x_values = [self.F(x) for x in range(len(values))]
        y_values = [self.F(y) for y in values]

        pol = [self.F(0)]
        k = len(x_values)
        for j in range(k):
            lj = [y_values[j]]
            for m in range(k):
                if m == j:
                    continue
                ljm = super()._mulPoly([1, - x_values[m]], [super()._RECIPROCAL(int(x_values[j] - x_values[m]))])
                lj = super()._mulPoly(lj, ljm)
            pol = self._addPoly(pol, lj)
        return _swap(pol)
    
    def eval_poly_at(self, coefficients: List[Field], index: Field):
        result = self.F(0)
        power_of_x = self.F(1)

        for coeff in coefficients:
            result = (result + (power_of_x * coeff))
            power_of_x = (power_of_x * index)

        return result
    
    def generate_commitment(self, coefficients):
        return super().generate_commitment(coefficients)
    
    def generate_proof(self, coefficients: List[Field], index: Field):
        quotientCoefficients = self._generate_quotient_polynomial(coefficients, index)
        #print("here", quotientCoefficients)
        return super().generate_commitment(quotientCoefficients)
    
    def _generate_quotient_polynomial(self, coefficients: List[Field], xVal: Field):
        yVal = self.eval_poly_at(coefficients, xVal)
        #print("Y-val quot:", yVal)
        x = [self.F(0), self.F(1)]
        res = super()._divPoly(super()._subPoly(coefficients, [yVal]), super()._subPoly(x, [xVal]))[0]
        #print("Res", res)
        return res

class Monomial(KZG10):
    def __init__(self, field=cu.curve_order) -> None:
        self.field = field
        self.F = GF(field)
    
    def interpolate(self, values: List[int]):
        x_values = [self.F(x) for x in range(len(values))]
        y_values = [self.F(y) for y in values]

        n = len(x_values)

        # Create the Vandermonde matrix
        A = []
        for i in range(n):
            row = [x_values[i]**j for j in range(n)]
            A.append(row)

        # Solve the system of equations to obtain coefficients
        coefficients = self.solve_system(A, y_values)

        return coefficients
    
    def solve_system(self, A, b):
        n = len(A)
        
        # Forward elimination (Gaussian elimination)
        for i in range(n):
            pivot = A[i][i]
            for j in range(i + 1, n):
                factor = A[j][i] / pivot
                for k in range(i, n):
                    A[j][k] -= factor * A[i][k]
                b[j] -= factor * b[i]

        # Backward substitution
        x = [0] * n
        for i in range(n - 1, -1, -1):
            x[i] = b[i]
            for j in range(i + 1, n):
                x[i] -= A[i][j] * x[j]
            x[i] /= A[i][i]

        return x
    
    def generate_commitment(self, coefficients):
        return super().generate_commitment(coefficients)
    

    def eval_monomial_at(self, coefficients, index):
        result = self.F(0)
        power_of_x = self.F(1)

        for coeff in coefficients:
            result += power_of_x * coeff
            power_of_x *= index

        return result
    
    def generate_proof(self, coefficients, index):
        quotient_coefficients = self._generate_quotient_polynomial(coefficients, index)
        return self.generate_commitment(quotient_coefficients)

    def _generate_quotient_polynomial(self, coefficients, x_val):
        y_val = self.eval_monomial_at(coefficients, x_val)
        x = [self.F(0), self.F(1)]
        res = super()._divPoly(super()._subPoly(coefficients, [y_val]), super()._subPoly(x, [x_val]))[0]
        return res

if __name__ == "__main__":
    main()
