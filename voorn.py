#!/usr/bin/env python

import random
from math import *
from numpy import *
from numpy.linalg import inv
import matplotlib.pyplot as plt
import StringIO
import mpld3
from mpld3 import plugins
from scipy.optimize import fsolve

from flask import Flask, request, make_response, render_template
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import json

""" Voorn-Overbeek """

def vorn_Spinodal(alpha,N):
		x = arange(1e-3,0.1,0.0001)
 		spinodal = ((2 * (2**.333) * ((N*x -x +1)**.666))/((3**.666)*(alpha**.666)*(N**.666)*(((x-1)**2)**(1./3.))*(x**.333)))
		return x, spinodal

def vjac(x,phi1,sigma,alpha,m):
	"df1/dphi2, df1/dchi; df2/dphi2, df2/dchi"
	"""
		1.0 - 2*m - 1.0/x[0] + m*1.0/(1 - x[0] - x[1]) 
		- (m*x[0])/(1.0 - x[0] - x[1]) - (m*x[1])/(1.0 - x[0] - x[1]) 
		+ (3.0*alpha*m*(sigma**2))/(4.0*(sigma*x[0] + x[1])**.5) 
		- (alpha*m*(sigma**2)*x[0])/(4*(sigma*x[0] + x[1])**0.5) 
		- (alpha*m*sigma*x[1])/(4.0*(sigma*x[0] + x[1])**.5) 
		- (0.5)*alpha*m*sigma*(sigma*x[0] + x[1])**.5, # dF1/dphi2

		-1.0*(m*1.0/(1 - phi1 - x[1])) + (m*phi1)/(1 - phi1 - x[1]) 
		+ m*1.0/(1 - x[0] - x[1]) - (m*x[0])/(1 - x[0] - x[1]) 
		+ (m*x[1])/(1 - phi1 - x[1]) - (m*x[1])/(1 - x[0] - x[1]) 
		- (3*alpha*m*sigma)/(4*(phi1*sigma + x[1])**.5) 
		+ (alpha*m*phi1*sigma)/(4*(phi1*sigma + x[1])**.5) 
		+ (alpha*m*x[1])/(4*(phi1*sigma + x[1])**.5) 
		+ (1.0/2.0)*alpha*m*(phi1*sigma + x[1])**.5 
		+ (3*alpha*m*sigma)/(4*(sigma*x[0] + x[1])**.5) 
		- (alpha*m*sigma*x[0])/(4*(sigma*x[0] + x[1])**.5) 
		- (alpha*m*x[1])/(4*(sigma*x[0] + x[1])**.5) 
		- (1.0/2.0)*alpha*m*(sigma*x[0] + x[1])**.5	
	"""
	return array([[
		 - (1./(m*x[0])) - 1./(1. - x[0]- x[1]) 
		+ (3*alpha*sigma**2)/(4*sqrt(x[1]+ x[0]*sigma)),

		
		1.0/(1 - phi1 - x[1]) - 1.0/(1 - x[0] - x[1]) 
		- (3*alpha*sigma)/(4.*sqrt(x[1] + phi1*sigma)) 
		+ (3*alpha*sigma)/(4.*sqrt(x[1] + x[0]*sigma))
		

		
		], #dF1/dpsi

		[log(phi1/2.0)/m - log(x[0]/2.0)/(m*1.0) - log(1 - phi1 - x[1]) 
		+ log(1 - x[0] - x[1]) - (1.5)*alpha*sigma*(phi1*sigma + x[1])**.5 
		+ (1.5)*alpha*sigma*(sigma*x[0] + x[1])**.5 , #dF2/dphi2

  		-log(1 - phi1 - x[1]) + log(1 - x[0] - x[1]) - (1.5)*alpha*(phi1*sigma 
		+ x[1])**.5 + (1.5)*alpha*(sigma*x[0] + x[1])**.5 
		+ (-phi1 + x[0])*(1.0/(1 - phi1 - x[1]) 
		- (3*alpha*sigma)/(4*(phi1*sigma + x[1])**.5))   ]]) #dF2/dpsi


def vfun(x,phi1,sigma,alpha,m):
	"F1 = f'(phi_1a) - f'(phi_2a); F2 = (b-a)*f'(phi_1a) -[ f(phi_2a) - f(phi_1a) ]"
	"""
		- phi1 + x[0] - m*(1 - phi1 - x[1]) + m*(1 - x[0] - x[1]) 
		- (1.5)*alpha*m*sigma*(phi1*sigma + x[1])**.5 
		+ (0.5)*alpha*m*phi1*sigma*(phi1*sigma + x[1])**0.5 
		+ 0.5*alpha*m*x[1]*(phi1*sigma + x[1])**0.5 
		+ (1.5)*alpha*m*sigma*(sigma*x[0] + x[1])**.5 
		- (0.5)*alpha*m*sigma*x[0]*(sigma*x[0] + x[1])**.5 
		- (0.5)*alpha*m*x[1]*(sigma*x[0] + x[1])**.5 + log(phi1/2.0) 
		- log(x[0]/2.0) + m*log(1 - phi1 - x[1]) - m*phi1*log(1 - phi1 - x[1]) 
		- m*(1 - phi1 - x[1])*log(1 - phi1 - x[1]) - m*x[1]*log(1 - phi1 - x[1]) 
		- m*log(1 - x[0] - x[1]) + m*x[0]*log(1 - x[0] - x[1]) 
		+ m*(1 - x[0] - x[1])*log(1 - x[0] - x[1]) + m*x[1]*log(1 - x[0] - x[1])

	"""
	return array([
		   (-1.5)*alpha*sigma*sqrt(x[1] + phi1*sigma) 
		+ (1.5)*alpha*sigma*sqrt(x[1] + x[0]*sigma) 
		+ log(phi1/2.)/m - log(x[0]/2.)/m - log(1 - phi1 - x[1]) 
		+ log(1 - x[0] - x[1]) 

		
		,

  		(-alpha)*(x[1] + phi1*sigma)**1.5 + alpha*(x[1] + x[0]*sigma)**1.5 
		+ (phi1*log(phi1/2))/m - (x[0]*log(x[0]/2))/m
		+ (-phi1 + x[0]) * (-1 + 1.0/m - 1.5*alpha*sigma*(x[1] + phi1*sigma)**.5 
		+ log(phi1/2)/m - log(1-phi1-x[1])) + (1-phi1-x[1])*log(1-phi1-x[1]) 
		- (1-x[0]-x[1])*log(1.-x[0]-x[1])
  ])

def v_crit(alpha,N):
		crit_phi = (-(N+2) + sqrt((N+2)**2 + 4*(N-1)))/(2*(N-1))
		crit_phi = crit_phi - .0001
		return crit_phi


def vNR(alpha,N,sigma):
		" Newton Raphson solver for the binary mixture"
		# Set up parameters, initial guesses, formatting, initializing etc.
		crit_phi = v_crit(alpha,N)
		print crit_phi

		phi1vals = arange(1e-5,.1,.002)
		phi1vals = phi1vals.tolist()
		guess = [0,0]
		new_guess = [0.1,.1] #phi2, psi
		iter = 0
		y2 = zeros((len(phi1vals),1))
		x2 = zeros((len(phi1vals),1))
		x1 = zeros((len(phi1vals),1))
		max_iter = 2000

		#Loop to find the roots using Multivariate Newton-Rhapson
		for phi in phi1vals:
			iter = 0
			while iter < max_iter :
				iter += 1
				index = phi1vals.index(phi)
				guess = new_guess
				jacobian = vjac(guess,phi,sigma,alpha,N)
				invjac = inv(jacobian)
				f1 = vfun(guess,phi,sigma,alpha,N)
				new_guess = guess - .1*dot(invjac,f1)
				if abs(new_guess[0] - guess[0]) < 1e-8 and abs(new_guess[1]-guess[1]) < 1e-8: 
					x1[index] = phi
					x2[index] = new_guess[0]
					y2[index] = new_guess[1]
					print x1[index], x2[index], y2[index]
					break

	#Convert Numpy arrays (x1,x2,y2) to a list
		x1=x1.tolist()
		x2=x2.tolist()
		x2=x2[::-1] #Has to reverse the order of x2, which was converted to a tuple in the previous line
		y2=y2.tolist()
		y2i = y2[::-1]

		#Concatenate the lists together
		phi = x1
		phi2 = x2

#		phi = x1 + x2
#		y2 = y2 + y2i
		return (phi,phi2,y2,y2i)

na = 1
nb = 1
crit_chi = .5 
crit_phi = 1
alpha = 3.655
N = 1000
sigma = .44

phi,phi2,y2,y2i = vNR(alpha,N,sigma)
plt.plot(phi,y2)
plt.plot(phi2,y2i)
plt.xlabel('phi')
plt.ylabel('psi')
plt.show()
#print phi,y2

# [phi2,psi]
x = [0.5,0.2]
phi1 = 0.2
m = 100.0



