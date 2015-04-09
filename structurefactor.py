#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import *
import numpy as np
import matplotlib.pyplot as plt

def structure_factor(na,nb,ba,bb,phi,chi):
	#Radius of Gyration
	Rga = na*ba**2 / 6.0 
	Rgb = nb*bb**2 / 6.0

	maxb = min(ba,bb)

	endq = 1.0/maxb

	qvals = np.arange(0.01,endq,0.01)
	sinv = np.zeros(( len(qvals) ))

	i = 0
	for q in qvals:

		#G_D(q^2,Rg^2)
		xa = q**2 * Rga**2
		g_a =(2.0/(xa**2)) * (xa - 1 + exp(-1*xa) )

		xb = q**2 * Rgb**2
		g_b =(2.0/(xb**2)) * (xb - 1 + exp(-1*xb) )

		#Inverse Structure Factor
		S1 = 1.0 / (na * phi * g_a)
		S2 = 1.0 / (nb * (1-phi) * g_b)
		sinv[i] = S1 + S2 - 2*chi
		sinv[i] = 1.0/sinv[i]
		i += 1
	return qvals,sinv

na = 1000
nb = 1000
ba = 2.2
bb = 2.1
phi = 0.8
chi = 0.8
"""
x = np.arange(0.0,5.0,0.1)
y = x**2
plt.plot(x,y)
plt.plot(y,x)
plt.show()
q,sinv = structure_factor(na,nb,ba,bb,phi,chi)
plt.plot(q,sinv)
plt.show()
"""

