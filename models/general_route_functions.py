import matplotlib.pyplot as plt
from math import *
import numpy as np
import SLCT

def flip(a1,a2,b1,b2,c1,c2):
		""" Switch a1 with a2, b1 with b2, c1 with c2"""
		return a2,a1,b2,b1,c2,c1

def vcrit():
	#Not sure where this is referenced, fyi,possible artifact from edits past
	firstderiv = (1 + log(phi))/N + (-1)*(1 + log(1-phi-psi)) - (3*sigma*alpha/2.)*(sigma*phi + psi)**(1./2.)
	secondderiv = 1./(N*phi) + 1./(1-phi-psi) - (3*sigma*sigma*alpha/4.)*(sigma*phi + psi)**(-1./2.)
	thirdderiv = -1./(N*phi*phi) + (1./(1-phi-psi))**2 + (3*alpha*(sigma**3)/8.)*(sigma*phi + psi)**(-3./2.)
