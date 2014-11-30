import numpy as np
from scipy import optimize

def fun(x):
		a = 1
		b = 1
		c = 1
		k = 1
		m = 100
		f = (-2*a*c*x**2)*(1-x)**2 + (2*c*(k-1)*x**3 + (b*(k-1) - c*(4*k - 1))*x**2 + 2*(c - b)*k*x + b*k)/(m*k)
		print f
		return f

a = 0.0
b = 1.3
c = 0.0
k = 1.0
m = 100
coeff = [2*a*c, 2*c*(k-1)/(m*k), (b*(k-1) - c*(4*k - 1))/(m*k) - 2*a*c , 2*(c - b)/m , b/m]

phi_c_temp =  np.roots(coeff)

for critval in phi_c_temp:
		if critval > 0 and critval < 1:
				phi_c = critval
				print phi_c

