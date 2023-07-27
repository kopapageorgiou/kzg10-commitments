from KZG10 import *
import numpy as np
import galois
def divided_diff(kzg: KZG10, x, y):
    n = len(y)
    coef = [[kzg.F(0)] * n for _ in range(n)]
    for i in range(n):
        coef[i][0] = y[i]
    
    for j in range(1, n):
        for i in range(n - j):
            coef[i][j] = (coef[i + 1][j - 1] - coef[i][j - 1]) / (x[i + j] - x[i])

    return coef


kzg = KZG10(43)

newton = kzg.choose_method(KZG10.NEWTON)

x1 = [0, 1, 2, 3]
y1 = [1, 3, 5, 9]
x2 = (0, 1, 2, 3)
y2 = (1, 3, 5, 9)
x = [kzg.F(i) for i in x1]
y = [kzg.F(i) for i in y1]
x1 = np.array(x)
y1 = np.array(y)
coeffs = newton.interpolate(y)
print(coeffs)
proof = newton.generate_proof(coeffs, 1)
coeffs = [29, 0, 2, 1]
GF = galois.GF(43)
f = galois.Poly([int(d) for d in coeffs], field=GF)
print("f=", f)
yPoly = galois.Poly([3], field=GF)
fmy = f - yPoly
print("fmy=", fmy)
xPoly = galois.Poly([1, 0], field=GF)
print("x=", xPoly)
xvalPoly = galois.Poly([1], galois.GF(43))
xmxval = xPoly - xvalPoly
print("xmxval=", xmxval)
#g = galois.Poly([-1, 1], field=GF)
divFG = fmy // xmxval
coeffs = divFG.coeffs
print(divFG)
print(coeffs)



