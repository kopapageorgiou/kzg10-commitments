from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

y = [5, 7, 12, 6]
x = [0, 1, 2, 3]
tck = interpolate.splrep(x, y, s=0, k=3) 
y_fit = interpolate.BSpline(*tck)(x)
print(dir(y_fit))
