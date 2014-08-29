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


app = Flask(__name__)


@app.route('/')

@app.route('/index')
def index():
    return render_template("index.html")

@app.route('/Flory.html',methods=['POST','GET'])
def flory():
	return render_template("Flory.html")

	
@app.route('/login', methods=['GET','POST'])	
def login():
	if request.method == 'POST':
		na = float(request.form['NFA'])
		nb = float(request.form['NFB'])

		crit_chi = .5*((1/(na**.5) + 1/(nb**.5))**2)
		nav = 2./crit_chi
 
		plt.xkcd()	
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		x = arange(0.1,0.90,0.001)
		spinodal = nav*(.5*(1./(na*x) + 1./(nb-nb*x)))
		phi,y2 =  NR(na,nb,nav)
		spinline = axis.plot(x,spinodal,'r',lw=2) 
		binline = axis.plot(phi,y2,'b',lw=2)
		fig.suptitle('Phase Diagram')
		canvas = FigureCanvas(fig)
		output = StringIO.StringIO()
		canvas.print_png(output, bbox_inches='tight')
		response = make_response(output.getvalue())
		response.mimetype = 'image/png'
		plugins.connect(fig, plugins.MousePosition())
		return mpld3.fig_to_html(fig)



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


def NR(na,nb,nav):
		" Newton Raphson solver for the binary mixture"
		# Set up parameters, initial guesses, formatting, initializing etc.

		if na != nb:
			crit_phi = (-nb + sqrt(na*nb))/(na-nb)
		else:
			crit_phi = .5  	

		phi1vals = arange(.001,crit_phi,.01)
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
		x1=x1.tolist()

		x2=x2.tolist()
		x2=x2[::-1] #Has to reverse the order of x2, which was converted to a tuple in the previous line
		y2 = nav*y2
		y2=y2.tolist()
		y2i = y2[::-1]

		#Concatenate the lists together
		phi = x1 + x2
		y2 = y2 + y2i
		return (phi,y2)
		#####################PLOT#######################
	

na = 1
nb = 1
crit_chi = .5 
crit_phi = 1

if __name__ == '__main__':
    app.run(debug=True)
#	app.run(host='0.0.0.0')
















