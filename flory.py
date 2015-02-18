#!/usr/bin/env python

import random
from math import *
from numpy import *
from numpy.linalg import inv
import matplotlib.pyplot as plt
import StringIO
import mpld3
from mpld3 import plugins

from flask import Flask, request, make_response, render_template
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import json

app = Flask(__name__)

@app.route('/')

@app.route('/index')
def index():
    return render_template("index.html")

@app.route('/howto.html')
def howto():
	return render_template("howto.html")

@app.route('/flory.html',methods=['POST','GET'])
def flory():
	return render_template("flory.html")

@app.route('/vorn.html',methods=['POST','GET'])
def vorn():
	return render_template("vorn.html")

@app.route('/slct.html',methods=['POST','GET'])
def slct():
	return render_template("slct.html")

@app.route('/slctplot', methods=['GET','POST'])	
def slctplot():
	if request.method == 'POST':
		na = float(request.form['NFA'])
		nb = float(request.form['NFB'])
		polya = request.form['polya']
		polyb = request.form['polyb']
		k1 = float(request.form['k1'])
		k2 = float(request.form['k2'])
		m1 = float(request.form['m1'])
		m2 = float(request.form['m2'])
		print polya,polyb,k1,k2,m1,m2

		z = 6.0
		
		""" Parameters for specific polymers"""
		#FIXMELATER
		""" Should encapsulate this elsewhere eventually """

		if polya == "PA":
				r1 = 1.0
				p1 = 1.0
		if polya == "PB":
				r1 = 1.333
				p1 = 1.333
		if polya == "PC":
				r1 = 6./5.
				p1 = 6./5.
		if polya == "PD":
				r1 = 4./3.
				p1 = 9./6.
		if polya == "PE":
				r1 = 7./4.
				p1 = 6./4.
		if polya == "PF":
				r1 = (3.+k1)/(2.+k1)
				p1 = (4.+k1)/(2.+k1)
		if polya == "PG":
				r1 = (5.+k1)/(4.+k1)
				p1 = (6.+k1)/(4.+k1)
		if polya == "PH":
				r1 = (4.+k1)/(3.+k1)
				p1 = (5.+k1)/(3.+k1)
		if polya == "PI":
				r1 = (6.+k1)/(3.+k1)
				p1 = (7.+k1)/(3.+k1)
		if polya == "PJ":
				r1 = (5.+k1+m1)/(2.+k1+m1)
				p1 = (8.+k1+m1)/(2.+k1+m1)
		if polya == "PK":
				r1 = 7./5.
				p1 = 8./5.
		if polya == "PL":
				r1 = 5./3.
				p1 = 5./3.
		if polya == "PM":
				r1 = 16./9.
				p1 = 19./9.
		if polya == "PN":
				r1 = 4./3.
				p1 = 5./3.
		if polya == "PO":
				r1 = 10./7.
				p1 = 12./7.
		if polya == "PP":
				r1 = 9./7.
				p1 = 14./7.
		if polya == "PQ":
				r1 = 11./7.
				p1 = 14./7.
		if polya == "PR":
				r1 = 11./7.
				p1 = 13./7.
		if polya == "PS":
				r1 = 13./8.
				p1 = 16./9.
		if polya == "PT":
				r1 = 11./8.
				p1 = 14./8.
		if polya == "PU":
				r1 = 13./9.
				p1 = 16./9.

		if polyb == "PA":
				r2 = 1.0
				p2 = 1.0
		if polyb == "PB":
				r2 = 1.333
				p2 = 1.333
		if polyb == "PC":
				r2 = 6./5.
				p2 = 6./5.
		if polyb == "PD":
				r2 = 4./3.
				p2 = 9./6.
		if polyb == "PE":
				r2 = 7./4.
				p2 = 6./4.
		if polyb == "PF":
				r2 = (3.+k2)/(2.+k2)
				p2 = (4.+k2)/(2.+k2)
		if polyb == "PG":
				r2 = (5.+k2)/(4.+k2)
				p2 = (6.+k2)/(4.+k2)
		if polyb == "PH":
				r2 = (4.+k2)/(3.+k2)
				p2 = (5.+k2)/(3.+k2)
		if polyb == "PI":
				r2 = (6.+k2)/(3.+k2)
				p2 = (7.+k2)/(3.+k2)
		if polyb == "PJ":
				r2 = (5.+k2+m2)/(2.+k2+m2)
				p2 = (8.+k2+m2)/(2.+k2+m2)
		if polyb == "PK":
				r2 = 7./5.
				p2 = 8./5.
		if polyb == "PL":
				r2 = 5./3.
				p2 = 5./3.
		if polyb == "PM":
				r2 = 16./9.
				p2 = 19./9.
		if polyb == "PN":
				r2 = 4./3.
				p2 = 5./3.
		if polyb == "PO":
				r2 = 10./7.
				p2 = 12./7.
		if polyb == "PP":
				r2 = 9./7.
				p2 = 14./7.
		if polyb == "PQ":
				r2 = 11./7.
				p2 = 14./7.
		if polyb == "PR":
				r2 = 11./7.
				p2 = 13./7.
		if polyb == "PS":
				r2 = 13./8.
				p2 = 16./9.
		if polyb == "PT":
				r2 = 11./8.
				p2 = 14./8.
		if polyb == "PU":
				r2 = 13./9.
				p2 = 16./9.
		if polyb == "PS":
				r2 = 13.0/8.0
				p2 = 16.0/9.0
		print r1, r2, z, p1, p2, na, nb
		global flipper

		"""Nice trick to fix asymmetry issue"""
		if na > nb:
				print "we gotta flip!"
				flipper = 1
				na, nb, r1, r2, p1, p2 =  flip(na,nb,r1,r2,p1,p2)
				print r1, r2, z, p1, p2, na, nb
		else:
				flipper = 0

		"""Set up the plot"""
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Volume Fraction')
		axis.set_ylabel('eps/kbT')
		axis.set_title('SLCT Phase Diagram')

		"""Run Optimization"""
		phi, y2 = SLCT_NR(r1,r2,z,p1,p2,na,nb)
		spinx,spiny = SLCT_Spinodal(r1,r2,z,p1,p2,na,nb)


		"""Plot"""
		spinline = axis.plot(spinx,spiny,'r',lw=2) 
		binline = axis.plot(phi,y2,'b',lw=2)

		"""Add d3 stuff"""
		canvas = FigureCanvas(fig)
		output = StringIO.StringIO()
		canvas.print_png(output, bbox_inches='tight')
		plugins.connect(fig, plugins.MousePosition())

		return mpld3.fig_to_html(fig,template_type='simple')

	
@app.route('/plot', methods=['GET','POST'])	
def plot():
	
	if request.method == 'POST':
		na = float(request.form['NFA'])
		nb = float(request.form['NFB'])

		"""Spinodal"""
		crit_chi = .5*((1/(na**.5) + 1/(nb**.5))**2)
		nav = 2./crit_chi
		global flipper

		"""Flipper"""
		if na > nb:
				flipper = 1
				na, nb, w,x,y,z=  flip(na,nb,1,1,1,1)
		else:
				flipper = 0

 
		"""Set up the plot"""
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Volume Fraction')
		axis.set_ylabel('Chi')
		axis.set_title('Flory-Huggins Phase Diagram')

		"""Run Optimization"""
		"""Need to move these lines it's own function"""
		x = arange(0.05,0.95,0.001)
		spinodal = nav*(.5*(1./(na*x) + 1./(nb-nb*x)))

		if flipper == 1:
			x = 1 - x

		phi,y2 =  NR(na,nb,nav,crit_chi)
		spinline = axis.plot(x,spinodal,'r',lw=2) 
		binline = axis.plot(phi,y2,'b',lw=2)

		"""Add d3 stuff"""
		canvas = FigureCanvas(fig)
		output = StringIO.StringIO()
		canvas.print_png(output, bbox_inches='tight')
		plugins.connect(fig, plugins.MousePosition())

		return mpld3.fig_to_html(fig,template_type='simple')

@app.route('/vornplot', methods=['GET','POST'])	
def vornplot():
	if request.method == 'POST':
		N = float(request.form['N'])
		sigma = float(request.form['sigma'])

		"""Set up the plot"""
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Volume Fraction, Phi')
		axis.set_ylabel('Salt Concentration, Psi')
		axis.set_title('Voorn-Overbeek Phase Diagram')
		canvas = FigureCanvas(fig)

		"""Move Spinodal Elsewhere"""
		phi,y2 =  vNR(alpha,N,sigma)
		#x, spinodal = vorn_Spinodal(alpha,N)
		#line1 = axis.plot(x,spinodal,'r',lw=2)
		spinline = axis.plot(phi,y2,'b',lw=2) 

		"""Make this organized like the other stuff"""
		plugins.connect(fig, plugins.MousePosition())
		return mpld3.fig_to_html(fig)

		


def flip(a1,a2,b1,b2,c1,c2):
		""" Switch a1 with a2, b1 with b2, c1 with c2"""
		return a2,a1,b2,b1,c2,c1


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

		phi1vals = arange(1e-2,.1,.002)
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
		return (phi,y2)
		


""" Flory Huggins"""

def fun(x,na,nb,phi1):
	"F1 = f'(phi_1a) - f'(phi_2a); F2 = (b-a)*f'(phi_1a) -[ f(phi_2a) - f(phi_1a) ]"

	na = 1.0*na
	nb = 1.0*nb
	return array([ 

			log(phi1) - log(x[0]) + x[1]*x[0]*(1-x[0])*na
			- x[1]*phi1*(1-phi1)*na + x[0] - phi1 - x[1]*na*(1-x[0]) + x[1]*na*(1-phi1)
			+ (na/nb)*(1-x[0]) - (na/nb)*(1-phi1),
			
			(x[0] - phi1)*(1./na - 1./nb + x[1] - 2*x[1]*phi1
			- log(1-phi1)/nb + log(phi1)/na) - ((x[0]/na)*log(x[0])
			+ ((1-x[0])/nb)*log(1-x[0]) + x[1]*x[0]*(1-x[0]))
			+ ((phi1/na)*log(phi1) + ((1-phi1)/nb)*log(1-phi1) + x[1]*(phi1)*(1-phi1))
			])	




def jac(x,na,nb,phi1):
	"df1/dphi2, df1/dchi; df2/dphi2, df2/dchi"
	na = 1.0*na
	nb = 1.0*nb

	"Make me better looking somehow"
	return array([
			[1 + x[1]*na - na/nb + x[1]*na*(1-x[0]) - 1.0/x[0] - x[1]*na*x[0],
			na*(1-phi1) - na*(1-phi1)*phi1 - na*(1-x[0]) + na*(1-x[0])*(x[0])],
			[log(phi1)/na -log(x[0])/na - log(1-phi1)/nb + log(1-x[0])/nb
			- 2*x[1]*phi1 + 2*x[1]*x[0], (x[0] - phi1)**2]])


def NR(na,nb,nav,crit_chi):
		" Newton Raphson solver for the binary mixture"
		# Set up parameters, initial guesses, formatting, initializing etc.

		if na != nb:
			crit_phi = (-nb + sqrt(na*nb))/(na-nb)
		else:
			crit_phi = .5  	

		phi1vals = arange(.001,crit_phi-.001,.01)
		phi1vals = phi1vals.tolist()
		guess = [0,0]
		new_guess = [0.5,3]
		iter = 0

		x1 = zeros((len(phi1vals),1))
		x2 = zeros((len(phi1vals),1))
		y2 = zeros((len(phi1vals),1))
		max_iter = 2000

		#Loop to find the roots using Multivariate Newton-Rhapson
		for phi in phi1vals:
			iter = 0
			while iter < max_iter :
				iter += 1
				index = phi1vals.index(phi)
				guess = new_guess
				jacobian = jac(guess,na,nb,phi)
				invjac = inv(jacobian)
				f1 = fun(guess,na,nb,phi)
				new_guess = guess - .1*dot(invjac,f1)
				if abs(new_guess[0] - guess[0]) < 1e-8 and abs(new_guess[1]-guess[1]) < 1e-8: 
					x1[index] = phi
					x2[index] = new_guess[0]
					y2[index] = new_guess[1]
					break

		# Flips the function back
		if flipper ==1:
			x1 = 1 - x1
			x2 = 1 - x2

		n = size(x1) + 1
		x1 = reshape(append(x1,crit_phi),(n,1))
		x1=x1.tolist()
		x2=x2.tolist()

		#Has to reverse the order of x2, which was converted to a tuple in the previous line
		x2=x2[::-1] 

		#Adds crit chi to the end of y2
		y2 = reshape(append(y2,crit_chi),(n,1))
		y2 = nav*y2
		y2=y2.tolist()
		y2i = y2[::-1]
		y2i.pop(0)

		#Concatenate the lists together
		phi = x1 + x2
		y2 = y2 + y2i


		return (phi,y2)
		#####################PLOT#######################


""" Simple Lattice Cluster """
def SLCT_crit(r1,r2,z,p1,p2,na,nb):
		#Use numpy root to calculate phi_c
		a = (r1 - r2)**2 / z**2
		b =((z-2)/2 + (1/z)*(-2*p1 + p2)) #Technically this is b/(eps/kt) which is factored out
		c = (3/z)*(p1 - p2) #Technically c/(eps/kt) 
		m = na
		k = nb*1.0/na
		coeff = [2*a*c, 2*c*(k-1)/(m*k), (b*(k-1) - c*(4*k - 1))/(m*k) - 2*a*c , 2*(c - b)/m , b/m]

		phi_c_temp =  roots(coeff)


		"Make sure that you pick the root that is real, positive and bounded by 0 and 1"
		for critval in phi_c_temp:
			if critval > 0 and critval < 1 and critval.imag == 0:
				phi_c = critval.real

		"Calculate the critical temperature"
		Tc = 2*(b + c*phi_c)/(1.0/(m*phi_c) + 1.0/(m*k*(1-phi_c)) - 2*a)
		return phi_c, Tc



def SLCT_Spinodal(r1,r2,z,p1,p2,na,nb):
		phi = arange(0.01,.99,0.001)
		a = (r1 - r2)**2 / z**2
		b =((z-2)/2 + (1.0/z)*(-2*p1 + p2)) #Technically this is b/(eps/kt) which is factored out
		c = (3.0/z)*(p1 - p2) #Technically c / (eps/kt) 
		f = (.5*(1./(na*phi) + 1./(nb-nb*phi)))
		spin1 = (f - a) / (b + c*phi)

		if flipper == 1:
				phi = 1 - phi

		return phi,spin1

def SLCT_fun(x,phi1,r1,r2,z,p1,p2,na,nb):
		a = (r1 - r2)**2 / z**2
		m1 = na*1.0
		m2 = nb*1.0
		
		return array([
				a*m1*(1-phi1) - m1*(1-phi1)/m2 - phi1 - a*m1*(1-phi1)*phi1 - a*m1*(1-x[0]) 
				+ m1*(1-x[0])/m2 + x[0] + a*m1*(1-x[0])*x[0] - m1*(1-phi1)*x[1]
				+ m1*(1-phi1)*phi1*x[1] + m1*(1-x[0])*x[1] - m1*(1-x[0])*x[0]*x[1]
				+ phi1*(1-phi1)*(1-phi1)*m1*p1*x[1]/z - x[0]*(1-x[0])*(1-x[0])*m1*p1*x[1]/z
				- m1*(1-phi1)*(1-phi1)*p1*x[1]/z + m1*(1-phi1)*(1-phi1)*phi1*p1*x[1]/z 
				+ phi1*phi1*(1-phi1)*m1*p2*x[1]/z - x[0]*x[0]*(1-x[0])*m1*p2*x[1]/z
				- 2*m1*(1-phi1)*phi1*p2*x[1]/z + m1*(1-phi1)*phi1*phi1*p2*x[1]/z + m1*p1*(1-x[0])*(1-x[0])*x[1]/z  
				+ 2*m1*p2*(1-x[0])*x[0]*x[1]/z - m1*p1*(1-x[0])*(1-x[0])*x[0]*x[1]/z 
				- m1*p2*(1-x[0])*x[0]*x[0]*x[1]/z + 0.5*m1*(1-phi1)*x[1]*z - 0.5*m1*(1-phi1)*phi1*x[1]*z
				- 0.5*m1*(1-x[0])*x[1]*z + 0.5*m1*(1-x[0])*x[0]*x[1]*z + log(phi1) - log(x[0]),

				- phi1/m1 + phi1/m2 + a*phi1*phi1 + x[0]/m1 - x[0]/m2 - 2*a*phi1*x[0] + a*x[0]*x[0] 
				- phi1*phi1*x[1] + 2*phi1*x[0]*x[1] - x[0]*x[0]*x[1] - 2*phi1*phi1*p1*x[1]/z
				+ 2*phi1*phi1*phi1*p1*x[1]/z + phi1*phi1*p2*x[1]/z - 2*phi1*phi1*phi1*p2*x[1]/z 
				+ 4*phi1*p1*x[0]*x[1]/z - 3*phi1*phi1*p1*x[0]*x[1]/z - 2*phi1*p2*x[0]*x[1]/z
				+ 3*phi1*phi1*p2*x[0]*x[1]/z -2*p1*x[0]*x[0]*x[1]/z + p2*x[0]*x[0]*x[1]/z 
				+ p1*x[0]*x[0]*x[0]*x[1]/z - p2*x[0]*x[0]*x[0]*x[1]/z 
				+ 0.5*phi1*phi1*x[1]*z - phi1*x[0]*x[1]*z + 0.5*x[0]*x[0]*x[1]*z 
				+ log(1-phi1)/m2 - x[0]*log(1-phi1)/m2 + x[0]*log(phi1)/m1 
				- log(1-x[0])/m2 + x[0]*log(1-x[0])/m2 -x[0]*log(x[0])/m1
				]) 


def SLCT_jac(x,phi1,r1,r2,z,p1,p2,na,nb):
	a = (r1 - r2)**2 / z**2
	m1 = na*1.0
	m2 = nb*1.0
	
	return array([[
			1 + a*m1 - m1/m2 + a*m1*(1-x[0]) - 1.0/x[0] - a*m1*x[0] - m1*x[1] - m1*(1-x[0])*x[1] + m1*x[0]*x[1] 
			- 2*m1*p1*(1-x[0])*x[1]/z + 2*m1*p2*(1-x[0])*x[1]/z - 2*m1*p1*(1-x[0])*(1-x[0])*x[1]/z 
			- 2*m1*p2*x[0]*x[1]/z + 4*m1*p1*(1-x[0])*x[0]*x[1]/z - 4*m1*p2*(1-x[0])*x[0]*x[1]/z 
			+ 2*m1*p2*x[0]*x[0]*x[1]/z + m1*x[1]*z/2 + 0.5*m1*(1-x[0])*x[1]*z - 0.5*m1*x[0]*x[1]*z,

			2*m1*phi1 - m1*phi1*phi1 - 2*m1*x[0] + m1*x[0]*x[0] + 4*m1*phi1*p1/z - 5*m1*phi1*phi1*p1/z 
			+ 2*m1*phi1*phi1*phi1*p1/z - 2*m1*phi1*p2/z + 4*m1*phi1*phi1*p2/z - 2*m1*phi1*phi1*phi1*p2/z
			- 4*m1*p1*x[0]/z + 2*m1*p2*x[0]/z + 5*m1*p1*x[0]*x[0]/z - 4*m1*p2*x[0]*x[0]/z - 2*m1*p1*x[0]*x[0]*x[0]/z 
			+ 2*m1*p2*x[0]*x[0]*x[0]/z - m1*phi1*z + 0.5*m1*phi1*phi1*z + m1*x[0]*z - 0.5*m1*x[0]*x[0]*z],
			[ 

			- 1.0/m2 - 2*a*phi1 + 1.0/(m2*(1-x[0])) + 2*a*x[0] - x[0]/(m2*(1-x[0])) 
			+ 2*phi1*x[1] - 2*x[0]*x[1] + 4*phi1*p1*x[1]/z
			- 3*phi1*phi1*p1*x[1]/z - 2*phi1*p2*x[1]/z + 3*phi1*phi1*p2*x[1]/z 
			- 4*p1*x[0]*x[1]/z + 2*p2*x[0]*x[1]/z + 3*p1*x[0]*x[0]*x[1]/z
			- 3*p2*x[0]*x[0]*x[1]/z - phi1*x[1]*z + x[0]*x[1]*z - log(1-phi1)/m2 
			+ log(phi1)/m2 + log(1-x[0])/m2 - log(x[0])/m1,

			-phi1*phi1 + 2*phi1*x[0] - x[0]*x[0] - 2*phi1*phi1*p1/z 
			+ 2*phi1*phi1*phi1*p1/z + phi1*phi1*p2/z - 2*phi1*phi1*phi1*p2/z
			+ 4*phi1*p1*x[0]/z - 3*phi1*phi1*p1*x[0]/z - 2*phi1*p2*x[0]/z 
			+ 3*phi1*phi1*p2*x[0]/z - 2*p1*x[0]*x[0]/z + p2*x[0]*x[0]/z
			+ p1*x[0]*x[0]*x[0]/z - p2*x[0]*x[0]*x[0]/z + phi1*phi1*z/2 
			- phi1*x[0]*z + x[0]*x[0]*z/2


			]])


"""insert df2/dphi'' and df2/db here """


def SLCT_NR(r1,r2,z,p1,p2,na,nb):
		" Newton Raphson solver for the binary mixture"
		# Set up parameters, initial guesses, formatting, initializing etc.
		phi_c, Tc = SLCT_crit(r1,r2,z,p1,p2,na,nb)
	
		phi1vals = arange(.01,phi_c,.009)
		phi1vals = phi1vals.tolist()
		guess = [0,0]
		new_guess = [0.88,2]
		iter = 0
		length = len(phi1vals)
		y2 = zeros((length,1))
		x2 = zeros((length,1))
		x1 = zeros((length,1))
		max_iter = 5000
		#Loop to find the roots using Multivariate Newton-Rhapson
		for phi in phi1vals:
			iter = 0
			damp = 0.5
			while iter < max_iter :

				iter += 1
				index = phi1vals.index(phi)
				guess = new_guess
				jacobian = SLCT_jac(guess,phi,r1,r2,z,p1,p2,na,nb)
				invjac = inv(jacobian)
				invjac = inv(jacobian)
				f1 = SLCT_fun(guess,phi,r1,r2,z,p1,p2,na,nb)
				new_guess = guess - damp*dot(invjac,f1)
				if abs(new_guess[0] - guess[0]) < 1e-6 and abs(new_guess[1]-guess[1]) < 1e-6: 
					x1[index] = phi
					x2[index] = new_guess[0]
					y2[index] = new_guess[1]
					break
		#Convert Numpy arrays (x1,x2,y2) to a list
		if flipper ==1:
			x1 = 1 - x1
			x2 = 1 - x2

		x1=x1.tolist()
		x2=x2.tolist()
		y2=y2.tolist()
		x2 = x2[::-1] #Has to reverse the order of x2, which was converted to a tuple in the previous line
		y2i = y2[::-1]
		phi = x1 + x2
		y2 = y2 + y2i

		return (phi,y2)

def vcrit():
	firstderiv = (1 + log(phi))/N + (-1)*(1 + log(1-phi-psi)) - (3*sigma*alpha/2.)*(sigma*phi + psi)**(1./2.)
	secondderiv = 1./(N*phi) + 1./(1-phi-psi) - (3*sigma*sigma*alpha/4.)*(sigma*phi + psi)**(-1./2.)
	thirdderiv = -1./(N*phi*phi) + (1./(1-phi-psi))**2 + (3*alpha*(sigma**3)/8.)*(sigma*phi + psi)**(-3./2.)



na = 1
nb = 1
crit_chi = .5 
crit_phi = 1
alpha = 3.655
N = 100



if __name__ == '__main__':
    app.run(debug=True)








