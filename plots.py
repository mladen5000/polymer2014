#!/usr/bin/env python
from math import *
from numpy import *
import matplotlib.pyplot as plt

"Flory huggins formula"
def flory_G(phi,na,nb,chi):
	"""Plots free energy"""
	enthalpy = chi*phi*(1-phi)
	entropy = phi/na * log(phi) + (1.-phi)/nb * log(1-phi) 
	f = phi/na * log(phi) + (1.-phi)/nb * log(1-phi) + chi*phi*(1-phi)
	print "The value is",f
	return enthalpy,entropy,f

phi   = arange(0.0001,1.0,0.001)
na = 100.0
nb = 100.0
chi = 0.03

h = zeros(( len(phi) , 1 ))
s = zeros(( len(phi) , 1 ))
g = zeros(( len(phi) , 1 ))

i = 0

for current_phi in phi:
	h[i],s[i],g[i]= flory_G(current_phi,na,nb,chi)
	i += 1

plt.plot(phi,h)
plt.plot(phi,s)
plt.plot(phi,g)
plt.show()

"""
plt.plot(phi,y2)
plt.plot(phi2,y2i)
plt.xlabel('phi')
plt.ylabel('psi')
plt.show()
"""




