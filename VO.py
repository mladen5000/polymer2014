#!/usr/bin/env python

#Import necessary modules
from math import *
from numpy import *
from numpy.linalg import inv
from scipy.optimize import root,fsolve,newton

"""
THIS IS THE MODULE FOR VOORN-OVERBEEK FILES
ALL VOORN OVERBEEK FUNCTIONS, WITH THE EXCEPTION OF PLOTTING
AND FRONTEND RELATED THINGS ARE PUT HERE
"""

def vCRIT(x,sigma,alpha,m):
	"""2nd and 3rd derivative of free energy function wrt phi,
	Used to determine the critical point of the function"""

	#Initialize
	phi = x[0]
	psi = x[1]

	#Second and Third Derivative of Free Energy with respect to phi
	f2 = 1.0/(m*phi) - (3*alpha*sigma**2)/(4.*sqrt(sigma*phi + psi)) + 1.0/(1-phi-psi)
	f3 =  -1.0/(m*phi*phi) + 3*alpha*sigma**3/(8.*(sigma*phi + psi)**1.5) + 1.0/(1-phi-psi)**2

	return array([f2,f3])

def vspin(x,phi,sigma,alpha,m):
	""" Function, fed into solver, in order to determine the spinodal curve."""

	#Initialize
	psi = x[0]

	#Second derivative of Free energy with respect to phi
	f2 = 1.0/(m*phi) - (3*alpha*sigma**2)/(4.*sqrt(sigma*phi + psi)) + 1.0/(1-phi-psi)

	return array([f2])

def vSpinodal(sigma,alpha,m):
	"""The actual Spinodal generating function which calls vspin
	1. Feed some initial values into the function
	2. Call the function
	3. Phivals picks the range of phi
	4. Given your selected phi, sigma; guess a psi, and try to minimize f2, once f2 is zero, return that psi
	5. Feed back into flory function
	"""
	#Range of Phi
	phivals = arange(1e-2,0.10,0.001)

	i=0
	xvals = zeros((len(phivals)))
	for phi in phivals:

		#Guess Parameters
		x0 = zeros((1))

		#Guess psi
		x0.fill(0.01)

		#Call the solver, return a value of psi
		xvals[i] = newton(vspin, x0, args = (phi,sigma, alpha, m))
		i += 1
	return phivals,xvals



def vCriticalpoint(sigma,alpha,m):
	"""Generates the Critical Value"""

	#Guess Parameters, [phi,psi]
	x0 = zeros((2))

	#Guess phi and psi, this is very picky
	x0.fill(0.01)

	#Call the solver, returns [phi,psi]
	x = fsolve(vCRIT, x0, args = (sigma, alpha, m))

	return x


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
	df1dphi = - (1./(m*x[0])) - 1./(1. - x[0]- x[1]) + (3*alpha*sigma**2)/(4*sqrt(x[1]+ x[0]*sigma))

	df1dpsi = 1.0/(1 - phi1 - x[1]) - 1.0/(1 - x[0] - x[1]) - (3*alpha*sigma)/(4.*sqrt(x[1] + phi1*sigma)) + (3*alpha*sigma)/(4.*sqrt(x[1] + x[0]*sigma))

	df2dphi = log(phi1/2.0)/m - log(x[0]/2.0)/(m*1.0) - log(1 - phi1 - x[1]) + log(1 - x[0] - x[1]) - (1.5)*alpha*sigma*(phi1*sigma + x[1])**.5 + (1.5)*alpha*sigma*(sigma*x[0] + x[1])**.5 

	df2dpsi = -log(1 - phi1 - x[1]) + log(1 - x[0] - x[1]) - (1.5)*alpha*(phi1*sigma + x[1])**.5 + (1.5)*alpha*(sigma*x[0] + x[1])**.5 + (-phi1 + x[0])*(1.0/(1 - phi1 - x[1]) - (3*alpha*sigma)/(4*(phi1*sigma + x[1])**.5)) 

	return array([ [df1dphi,df1dpsi],[df2dphi,df2dpsi] ])


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
	

def vNR(alpha,N,sigma):
		""" Newton Raphson solver for the binary mixture"""
		# Set up parameters, initial guesses, formatting, initializing etc.

		critphi = vCriticalpoint(sigma,alpha,N)
		phi1vals = arange(1e-2,critphi[0],.002)
		phi1vals = phi1vals.tolist()
		guess = [0,0]
		new_guess = [0.1,0.1] #phi2, psi
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
				new_guess = guess - 0.1*dot(invjac,f1)
				if abs(new_guess[0] - guess[0]) < 1e-5 and abs(new_guess[1]-guess[1]) < 1e-5: 
					x1[index] = phi
					x2[index] = new_guess[0]
					y2[index] = new_guess[1]
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

		#This is the correct way to string the data together
		phi = x1 + x2
		psi = y2 + y2i

		return (phi,psi)

def vfree(N,psi,sigma):
	"""Calculates the free energy, enthalpy, entropy of VO"""
	#Define
	alpha = 3.655
	phivals = arange(0.0,1.0,0.001)

	#Initalize
	enthalpy = zeros(( len(phivals) ))
	entropy = zeros(( len(phivals) ))
	f = zeros(( len(phivals) ))
	
	#Loop and calculate values
	i=0
	for phi in phivals:
		enthalpy[i] = -alpha*(phi*sigma + psi)**1.5
		entropy[i] = (phi/N)*log(phi/2.0) + psi*log(psi/2.0) + (1-phi-psi)*log(1-phi-psi) 
		f[i] = (phi/N)*log(phi/2.0) + psi*log(psi/2.0) + (1-phi-psi)*log(1-phi-psi) - alpha*(phi*sigma + psi)**1.5
		i += 1

	return phivals, enthalpy,entropy, f



