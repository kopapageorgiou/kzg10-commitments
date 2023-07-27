from KZG10 import *

kzg = KZG10()
x = [0, 1, 2]
y = [1, 5, 17]

mon = kzg.choose_method(KZG10.MONOMIAL)
coeffs = mon.interpolate(y)
print(coeffs)