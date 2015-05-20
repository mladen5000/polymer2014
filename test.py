from math import *
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt


#I am currently using this to see if I can get a more advanced version of semiflexibility up and running
#Current one online assumes a constant temperature(obviously wrong)
#Trying to numerically solve for the Spinodal



def SLCT_flexspin(T,na,nb,flex1,flex2,eps,phi):
	"""Auxillary function, which is called by solve, finds the value of the Spinodal at a fixed phi"""

	z = 6.0

	
	#Parameters for A
	a0_a = flex1[0]
	a1_a = flex1[1]
	g111_a = flex1[2]
	g21_a = flex1[3]
	g3_a = flex1[4]
	Eb_a = flex1[5]

	#Parameters for B
	a0_b = flex2[0]
	a1_b = flex2[1]
	g111_b = flex2[2]
	g21_b = flex2[3]
	g3_b = flex2[4]
	Eb_b = flex2[5]


	#Dependent on T
	#A terms
	g_a = z /  (z - 1 + exp( Eb_a/T ) )
	r1 = a0_a + a1_a*g_a
	p1 = g111_a + g21_a*g_a + g3_a*(g_a**2)

	#B terms
	g_b = z /  (z - 1 + exp(Eb_b/T))
	r2 = a0_b + a1_b*g_b
	p2 = g111_b + g21_b*g_b + g3_b*g_b**2
	
	a = (r1 - r2)**2 / z**2
	b =((z-2)/2 + (1.0/z)*(-2*p1 + p2)) #Technically this is b/(eps/kt) which is factored out
	c = (3.0/z)*(p1 - p2) #Technically c / (eps/kt) 
	f = (.5*(1./(na*phi) + 1./(nb-nb*phi)))

	#Calculate the residual
	rhs = (f - a) / (b+c*phi)
	eps = 2.0
	res = eps/T - (f - a) / (b + c*phi)
	return res

def run_SLCT_flexspinodal(na,nb,flex1,flex2,eps,phi):
	phivals = np.arange(0.02,0.98,0.01)
	tempforphi = np.zeros(( len(phivals) ))

	i=0
	for phi in phivals:
		tempforphi[i] = fsolve(SLCT_flexspin,x0,args=(na,nb,flex1,flex2,eps,phi) ,factor=0.005)
		i = i+1
	return phivals, tempforphi

na = 100.0
nb = 100.0
flex1 = [0.0,1.0,0.0,0.0,1.0,1900.0]
flex2 = [0.0,1.0,0.0,0.0,1.0,00.0]
x0 = 485.0
eps = 2.0

i=0
tempforphi = np.zeros(( len(phivals) ))
for phi in phivals:
	tempforphi[i] = fsolve(SLCT_flexspin,x0,args=(na,nb,flex1,flex2,eps,phi) ,factor=0.005)
	i = i+1
	print phi, tempforphi

plt.plot(phivals,tempforphi)
plt.show()




