#!/usr/bin/env python

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
		req_json = request.get_json()
		model = req_json['model']
		na1 = req_json['NFA']
		na2 = req_json['NFB']
		na = float(na1)
		nb = float(na2)

		crit_chi = .5*((1/(na**.5) + 1/(nb**.5))**2)
		nav = 2./crit_chi
 
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Volume Fraction')
		axis.set_ylabel('Chi')
		axis.set_title('Flory-Huggins Phase Diagram')
		x = arange(0.05,0.95,0.001)
		axis = arange(0.00,1.00,0.001)
		spinodal = nav*(.5*(1./(na*x) + 1./(nb-nb*x)))
		phi,y2 =  NR(na,nb,nav,crit_chi)  # data
		
		data = list()
		header = list()
		header.append('x')
		header.append('spineline')
		header.append('y2')
		data.append(header)
		
		for a in range(0,len(phi)):
			if float(phi[a][0]) < 0.05:
				point = list()
				point.append(float(phi[a][0]))
				point.append("null")
				point.append(float(y2[a][0]))
				data.append(point)

			
		for i in range(0,len(x)):
			found_match = 0
			match = -1
			for j in range(0,len(phi)):
				test = ceil(phi[j][0] * 1000) / 1000.0
				test2 = ceil(x[i] * 1000) / 1000.0
				if test == test2:
					found_match = 1
					match = j
			
			if found_match != 0:
				point = list()
				point.append(float(str(x[i]).replace("[","").replace("]","")))
				test0 = ceil(spinodal[i] * 1000) / 1000.0
				point.append(float(str(test0).replace("[","").replace("]","")))
				test3 = ceil(y2[match][0] * 1000) / 1000.0
				point.append(float(str(test3).replace("[","").replace("]","")))
				data.append(point)
				
			else:
				point = list()
				point.append(float(str(x[i]).replace("[","").replace("]","")))
				test0 = ceil(spinodal[i] * 1000) / 1000.0
				point.append(float(str(test0).replace("[","").replace("]","")))
				point.append("null")
				data.append(point)
		
		for a in range(0,len(phi)):
			if float(phi[a][0]) > 0.95:
				point = list()
				point.append(float(phi[a][0]))
				point.append("null")
				point.append(float(y2[a][0]))
				data.append(point)
		
		return str(data).replace('"null"','null')

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
	

na = 1
nb = 1
crit_chi = .5 
crit_phi = 1
alpha = 3.655
N = 1





if __name__ == '__main__':
    app.run(debug=True)
















