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

def generate_SLCTfigure(na,nb,polya,polyb,k1,k2,m1,m2,eps):
		fig = plt.figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Volume Fraction, \u03a6')
		axis.set_ylabel('Free Energy')
		axis.set_title('SLCT Free Energy Diagram')

		"""Run Optimization"""
		"""Need to move these lines it's own function"""
		r1, p1, r2, p2 = SLCT.SLCT_constants(polya,polyb,k1,k2,m1,m2)
		z = 6.0
		phi,h,s,g = SLCT.SLCTfree(r1,r2,z,p1,p2,na,nb,eps)

		
		hline = axis.plot(phi,h,'r',lw=1,alpha = 0.5,label='Enthalpy')
		sline = axis.plot(phi,s,'b',lw=1,alpha = 0.5,label='Entropy')
		gline = axis.plot(phi,g,'g',lw=3,label="Free Energy")
		legend = axis.legend()
		return fig
