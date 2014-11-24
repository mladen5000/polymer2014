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

@app.route('/Flory.html',methods=['POST','GET'])
def flory():
	return render_template("Flory.html")

@app.route('/vorn.html',methods=['POST','GET'])
def vorn():
	return render_template("vorn.html")

	
@app.route('/plot', methods=['GET','POST'])	
def login():
	if request.method == 'POST':
		na = float(request.form['NFA'])
		nb = float(request.form['NFB'])

		crit_chi = .5*((1/(na**.5) + 1/(nb**.5))**2)
		nav = 2./crit_chi
 
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Volume Fraction')
		axis.set_ylabel('Chi')
		axis.set_title('Flory-Huggins Phase Diagram')
		x = arange(0.05,0.95,0.001)
		spinodal = nav*(.5*(1./(na*x) + 1./(nb-nb*x)))

		phi,y2 =  NR(na,nb,nav,crit_chi)
		spinline = axis.plot(x,spinodal,'r',lw=2) 
		binline = axis.plot(phi,y2,'b',lw=2)
		canvas = FigureCanvas(fig)
		output = StringIO.StringIO()
		canvas.print_png(output, bbox_inches='tight')
		plugins.connect(fig, plugins.MousePosition())

#		json01 = json.dumps(mpld3.fig_to_dict(fig,template_type='simple'))
		return mpld3.fig_to_html(fig,template_type='simple')
#		return render_template("plot.html",json01=json01)

@app.route('/vornplot', methods=['GET','POST'])	
def vornplot():
	if request.method == 'POST':
		N = float(request.form['N'])
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		x = arange(1e-5,0.1,0.0001)
 		spinodal = ((2 * (2**.333) * ((N*x -x +1)**.666))/((3**.666)*(alpha**.666)*(N**.666)*(((x-1)**2)**(1./3.))*(x**.333)))
		phi,y2 =  vNR(alpha,N)
		line1 = axis.plot(x,spinodal,'r',lw=2)
		spinline = axis.plot(phi,y2,'b',lw=2) 
		axis.set_xlabel('Volume Fraction')
		axis.set_ylabel('Charge Density')
		axis.set_title('Voorn-Overbeek Phase Diagram')
		canvas = FigureCanvas(fig)
#		output = StringIO.StringIO()
#		canvas.print_png(output, bbox_inches='tight')
		plugins.connect(fig, plugins.MousePosition())
		return mpld3.fig_to_html(fig)






""" Voorn-Overbeek """

def vfun(x,alpha,N,phi1):
	"F1 = f'(phi_1a) - f'(phi_2a); F2 = (b-a)*f'(phi_1a) -[ f(phi_2a) - f(phi_1a) ]"
	return array([1.5*alpha*x[1]*(x[1]*phi1)**0.5 - 1.5*alpha*x[1]*(x[1]*x[0])**0.5
			- log(phi1/2.)/N + log(x[0]/2.)/N + log(1-phi1) - log(1-x[0]),
			-1.5*alpha*x[1]*x[0]*(x[1]*phi1)**.5 + .5*alpha*x[1]*phi1*(x[1]*phi1)**.5 + alpha*x[1]*x[0]*(x[1]*x[0])**.5 + x[0]*log(phi1/2)/N - phi1/N + x[0]/N - x[0]*log(x[0]/2)/N - x[0]*log(1-phi1) + phi1 + log(1-phi1) - x[0] + x[0]*log(1-x[0]) - log(1-x[0])])

def vjac(x,alpha,N,phi1):
	"df1/dphi2, df1/dchi; df2/dphi2, df2/dchi"
	return array([[((-3.*alpha*x[1]**2.)/(4.*(x[1]*x[0])**0.5)) + 1./(N*x[0]) + 1./(1.-x[0]), # dF1/dphi2
		 2.25*alpha*((x[1]*phi1)**0.5 - (x[1]*x[0])**0.5)], #dF1/dsigma
		 [-1.5*alpha*x[1]*(x[1]*phi1)**0.5 + 1.5*alpha*x[1]*(x[1]*x[0])**0.5 + log(phi1)/N - log(x[0])/N - log(1.-phi1) + log(1.-x[0]), #dF2/dphi2
			0.75*alpha*(-3.*x[0]*(x[1]*phi1)**.5 + phi1*(x[1]*phi1)**.5 + 2.*x[0]*(x[1]*x[0])**.5)]]) #dF2/dsigma

def vNR(alpha,N):
		" Newton Raphson solver for the binary mixture"
		# Set up parameters, initial guesses, formatting, initializing etc.
		crit_phi = (-(N+2) + sqrt((N+2)**2 + 4*(N-1)))/(2*(N-1))
		crit_phi = crit_phi - .0001
		phi1vals = arange(2e-7,crit_phi,.0001)
		phi1vals = phi1vals.tolist()
		guess = [0,0]
		new_guess = [0.9,.9] #phi2, sigma
		iter = 0
		length = len(phi1vals)
		y2 = zeros((length,1))
		x2 = zeros((length,1))
		x1 = zeros((length,1))
		max_iter = 2000

		#Loop to find the roots using Multivariate Newton-Rhapson
		for phi in phi1vals:
			iter = 0
			while iter < max_iter :
				iter += 1
				index = phi1vals.index(phi)
				guess = new_guess
				jacobian = vjac(guess,alpha,N,phi)
				invjac = inv(jacobian)
				f1 = vfun(guess,alpha,N,phi)
				new_guess = guess - .1*dot(invjac,f1)
				if abs(new_guess[0] - guess[0]) < 1e-8 and abs(new_guess[1]-guess[1]) < 1e-8: 
					x1[index] = phi
					x2[index] = new_guess[0]
					y2[index] = new_guess[1]
					break

	#Convert Numpy arrays (x1,x2,y2) to a list
		x1=x1.tolist()
		x2=x2.tolist()
		print x1
		print x2
		print y2
		x2=x2[::-1] #Has to reverse the order of x2, which was converted to a tuple in the previous line
		y2=y2.tolist()
		y2i = y2[::-1]

		#Concatenate the lists together
		phi = x1 + x2
		y2 = y2 + y2i
		return (phi,y2)
		


""" Flory Huggins"""

def fun(x,na,nb,phi1):
	"F1 = f'(phi_1a) - f'(phi_2a); F2 = (b-a)*f'(phi_1a) -[ f(phi_2a) - f(phi_1a) ]"
	return array([log(x[0])/na - log(phi1)/na + log(1.-phi1)/nb -
			log(1.-x[0])/nb + 2.*x[1]*phi1 - 2.*x[1]*x[0],
			(x[0] - phi1)*(1./na - 1./nb + x[1] - 2*x[1]*phi1
			- log(1-phi1)/nb + log(phi1)/na) - ((x[0]/na)*log(x[0])
			+ ((1-x[0])/nb)*log(1-x[0]) + x[1]*x[0]*(1-x[0]))
			 + ((phi1/na)*log(phi1) + ((1-phi1)/nb)*log(1-phi1) + x[1]*(phi1)*(1-phi1))])	

def jac(x,na,nb,phi1):
	"df1/dphi2, df1/dchi; df2/dphi2, df2/dchi"
	return array([[1./(na*x[0]) + 1./(nb*(1.-x[0])) - 2.*x[1],  2.*(phi1 - x[0])],
				[log(phi1)/na -log(x[0])/na - log(1-phi1)/nb + log(1-x[0])/nb
				-2*x[1]*phi1 + 2*x[1]*x[0],
				(x[0] - phi1)**2]])


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
		new_guess = [0.3,3]
		iter = 0
		length = len(phi1vals)
		y2 = zeros((length,1))
		x2 = zeros((length,1))
		x1 = zeros((length,1))
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
		#Convert Numpy arrays (x1,x2,y2) to a list
		x1 = reshape(append(x1,crit_phi),(51,1))
		x1=x1.tolist()
		x2=x2.tolist()

		x2=x2[::-1] #Has to reverse the order of x2, which was converted to a tuple in the previous line
		y2 = reshape(append(y2,crit_chi),(51,1))
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

		for critval in phi_c_temp:
			if critval > 0 and critval < 1:
				phi_c = critval
				return phi_c



def SLCT_Spinodal(r1,r2,z,p1,p2,na,nb):
		#put all these values into a list
		phi = arange(0.01,.99,0.001)
		a = (r1 - r2)**2 / z**2
		b =((z-2)/2 + (1/z)*(-2*p1 + p2)) #Technically this is b/(eps/kt) which is factored out
		c = (3/z)*(p1 - p2) #Technically c / (eps/kt) 
		f = (.5*(1./(na*phi) + 1./(nb-nb*phi)))
		#spin1 = (f - a) / (b + c*phi)
		btmspin = (z - 2 + (2./z)*(p1*(1 - 3*(1-phi)) + p2*(1-3*phi)))
		topspin = (1./(na*phi) + 1./(nb*(1-phi)) - 2*a)
		spin2 = topspin/btmspin
		return phi,spin2

def SLCT_fun(x,phi1,r1,r2,z,p1,p2,na,nb):
		a = (r1 - r2)**2 / z**2
		c = (z-2)/2 
		d = 1.0/z
		m1 = na 
		m2 = nb
		
		k = nb*1.0/na
		
		return array([
		1.0/m1 - 1.0/(k*m1) - 1.0/m2 + 1.0/(k*m2) - 2*a*phi1 - 2*x[1]*c*phi1
		+ 4*x[1]*d*p1*phi1 - 2*x[1]*d*p2*phi1 
		- 3*x[1]*d*p1*phi1**2 + 3*x[1]*d*p2*phi1**2 + 2*a*x[0] + 2*x[1]*c*x[0]
		- 4*x[1]*d*p1*x[0] + 2*x[1]*d*p2*x[0] + 3*x[1]*d*p1*x[0]**2 - 3*x[1]*d*p2*x[0]**2
		- log(1-phi1)*(1.0/(k*m1)) + log(phi1)*(1.0/m1) + log(1-x[0])*(1.0/(k*m2)) - log(x[0])*(1.0/m2),

		phi1/m1 - phi1/(k*m1) -a*phi1**2 - x[1]*c*phi1**2 + 2*x[1]*d*p1*phi1**2
		- x[1]*d*p2*phi1**2 - 2*x[1]*d*p1*phi1**3 +2*x[1]*d*p2*phi1**3 - x[0]/m1
		+ x[0]/(k*m1)  + 2*a*phi1*x[0] + 2*x[1]*c*phi1*x[0] - 4*x[1]*d*p1*phi1*x[0]
		+ 2*x[1]*d*p2*phi1*x[0] + 3*x[1]*d*p1*x[0]*phi1**2 - 3*x[1]*d*p2*x[0]*phi1**2
		- a*x[0]**2 - x[1]*c*x[0]**2 + 2*x[1]*d*p1*x[0]**2 - x[1]*d*p2*x[0]**2
		- x[1]*d*p1*x[0]**3 + x[1]*d*p2*x[0]**3 - log(1-phi1)/(k*m1)
		+ x[0]*log(1-phi1)*(1.0/(k*m1)) - x[0]*log(phi1)*(1.0/m1) + log(1-x[0])*(1.0/(k*m2))
		- x[0]*log(1-x[0])*(1.0/(k*m2)) + x[0]*log(x[0])*(1.0/m2)
		]) 

	
def SLCT_jac(x,phi1,r1,r2,z,p1,p2,na,nb):
	a = (r1 - r2)**2 / z**2
	c = (z-2.0)/2
	d = 1.0/z
	m1 = na 
	m2 = nb
	k = nb*1.0/na
	
	return array([[
			2*a + 2*x[1]*c - 4*x[1]*d*p1 + 2*x[1]*d*p2
			- 1.0/(k*m2*(1-x[0])) - 1.0/(m2*x[0]) + 6*x[1]*d*p1*x[0]
			- 6*x[1]*d*p2*x[0],

			-2*c*phi1 + 4*d*p1*phi1 - 2*d*p2*phi1 - 3*d*p1*phi1**2
			+ 3*d*p2*phi1**2 + 2*c*x[0] - 4*d*p1*x[0] + 2*d*p2*x[0] 
			+ 3*d*p1*x[0]**2 - 3*d*p2*x[0]**2],

			[-1.0/m1 + 1.0/(k*m1) + 1.0/m2 - 1.0/(k*m2) + 2*a*phi1
		    + 2*x[1]*c*phi1 - 4*x[1]*d*p1*phi1
			+ 2*x[1]*d*p2*phi1 + 3*x[1]*d*p1*phi1**2 - 3*x[1]*d*p2*phi1**2
			- 2*a*x[0] - 2*x[1]*c*x[0] + 4*x[1]*d*p1*x[0]
			- 2*x[1]*d*p2*x[0]  - 3*x[1]*d*p1*x[0]**2
			+ 3*x[1]*d*p2*x[0]**2 + log(1-phi1)/(k*m1) - log(phi1)/m1 
			- log(1-x[0])*(1.0/(k*m2)) + log(x[0])*(1.0/m2),

			-c*phi1**2 + 2*d*p1*phi1**2 - d*p2*phi1**2 - 2*d*p1*phi1**3
			+ 2*d*p2*phi1**3 + 2*c*phi1*x[0] - 4*d*p1*phi1*x[0] + 2*d*p2*phi1*x[0]
			+ 3*d*p1*x[0]*phi1**2 - 3*d*p2*x[0]*phi1**2 - c*x[0]**2
			+ 2*d*p1*x[0]**2 - d*p2*x[0]**2 - d*p1*x[0]**3 + d*p2*x[0]**3
			]])


def SLCT_NR(r1,r2,z,p1,p2,na,nb):
		" Newton Raphson solver for the binary mixture"
		# Set up parameters, initial guesses, formatting, initializing etc.
		phi_c = SLCT_crit(r1,r2,z,p1,p2,na,nb)
		print phi_c
		fil1=open('./blah.txt', 'w+')
	
		phi1vals = arange(.01,phi_c,.01)
		phi1vals = phi1vals.tolist()
		guess = [0,0]
		new_guess = [0.99,.01]
		iter = 0
		length = len(phi1vals)
		y2 = zeros((length,1))
		x2 = zeros((length,1))
		x1 = zeros((length,1))
		max_iter = 2000
		#Loop to find the roots using Multivariate Newton-Rhapson
		for phi in phi1vals:
			iter = 0
			damp = 0.1
			while iter < max_iter :

				print >> fil1, guess[0], guess[1] 
				iter += 1
				index = phi1vals.index(phi)
				guess = new_guess
				if guess[0] < 0 or guess[0] > 1:
						guess[0] = random.random()
						guess[1] = 0
						damp = 0.01
						#print phi, iter, damp, guess[0]
				jacobian = SLCT_jac(guess,phi,r1,r2,z,p1,p2,na,nb)
				invjac = inv(jacobian)
				invjac = inv(jacobian)
				f1 = SLCT_fun(guess,phi,r1,r2,z,p1,p2,na,nb)
				new_guess = guess - damp*dot(invjac,f1)
				if abs(new_guess[0] - guess[0]) < 1e-14 and abs(new_guess[1]-guess[1]) < 1e-14: 
					x1[index] = phi
					x2[index] = new_guess[0]
					y2[index] = new_guess[1]
					break
		#Convert Numpy arrays (x1,x2,y2) to a list
		print x1
		print x2
		print y2
		x1=x1.tolist()
		x2=x2.tolist()
		y2=y2.tolist()
		x2 = x2[::-1] #Has to reverse the order of x2, which was converted to a tuple in the previous line
		y2i = y2[::-1]
		phi = x1 + x2
		y2 = y2 + y2i
		return (phi,y2)





na = 1
nb = 1
crit_chi = .5 
crit_phi = 1
alpha = 3.655
N = 1
phi, y2 = SLCT_NR(  1.75,1.2,4,1.5,1.2,100,100)
phix,spinx2 = SLCT_Spinodal(1.75,1.2,4,1.5,1.2,100,100)
plt.plot(phi,y2)
plt.plot(phix,spinx2)
"""
x = [0.5,0.1]
phi1 = 0.1
print SLCT_jac(x,phi1,1.75,1.2,4,1.5,1.2,100,100)
"""

plt.show()
if __name__ == '__main__':
    app.run(debug=True)
















