#!/usr/bin/env python

import random
from math import *
from numpy.linalg import inv
import numpy as np
import matplotlib.pyplot as plt
import StringIO
import mpld3
from mpld3 import plugins
from SLCT import *

import json
from general_route_functions import generate_figure
from flask import request

""" Flory Huggins"""
def fun(x,na,nb,phi1):
	"""These are the functions necessary to solve the binodal:
	F1 is the criteria that the chemical potentials are equal, whereas F2 is the definition of the derivative w.r.t. phi
	F1 = f'(phi_1a) - f'(phi_2a);
	F2 = (b-a) * f'(phi_1a) - [ f(phi_2a) - f(phi_1a) ]
	"""

	#Convert to floats
	na = 1.0*na
	nb = 1.0*nb

	#Actual expression of F1, and F2 are housed in this array, this is more for formatting reasons
	return np.array([ 

			log(phi1) - log(x[0]) + x[1]*x[0]*(1-x[0])*na
			- x[1]*phi1*(1-phi1)*na + x[0] - phi1 - x[1]*na*(1-x[0]) + x[1]*na*(1-phi1)
			+ (na/nb)*(1-x[0]) - (na/nb)*(1-phi1),

			(x[0] - phi1)*(1./na - 1./nb + x[1] - 2*x[1]*phi1
			- log(1-phi1)/nb + log(phi1)/na) - ((x[0]/na)*log(x[0])
			+ ((1-x[0])/nb)*log(1-x[0]) + x[1]*x[0]*(1-x[0]))
			+ ((phi1/na)*log(phi1) + ((1-phi1)/nb)*log(1-phi1) + x[1]*(phi1)*(1-phi1))
			])	


def jac(x,na,nb,phi1):
	""" This the expression for the jacobian, given by the following:
	[[df1/dphi2, df1/dchi]
	[df2/dphi2, df2/dchi] ]
	"""
	#Convert to floats
	na = 1.0*na
	nb = 1.0*nb

	return np.array([
			[1 + x[1]*na - na/nb + x[1]*na*(1-x[0]) - 1.0/x[0] - x[1]*na*x[0],
			na*(1-phi1) - na*(1-phi1)*phi1 - na*(1-x[0]) + na*(1-x[0])*(x[0])],
			[log(phi1)/na -log(x[0])/na - log(1-phi1)/nb + log(1-x[0])/nb
			- 2*x[1]*phi1 + 2*x[1]*x[0], (x[0] - phi1)**2]])


def NR(na,nb,nav,crit_chi,flipper):
		" Newton Raphson solver for the binary mixture"
		# Set up parameters, initial guesses, formatting, initializing etc.

		#Establish the critical point analytically
		if na != nb:
			crit_phi = (-nb + sqrt(na*nb))/(na-nb)
		else:
			crit_phi = .5  	

		#Set up the array of phi_a in phase 1
		phi1vals = np.arange(.001,crit_phi-.001,.01)
		phi1vals = phi1vals.tolist()

		#Establish initial guess
		guess = [0,0]
		new_guess = [0.5,3]
		iter = 0

		#Generate Arrays
		x1 = np.zeros((len(phi1vals),1)) # Final array to hold phi_a in phase 1
		x2 = np.zeros((len(phi1vals),1)) # Final array to hold phi_a in phase 2
		y2 = np.zeros((len(phi1vals),1)) # Final array to hold chi

		max_iter = 2000
		damp = 0.1 #Damping constant to help the solver 

		#Loop to find the roots using Multivariate Newton-Rhapson
		for phi in phi1vals:
			iter = 0
			while iter < max_iter :
				iter += 1
				index = phi1vals.index(phi)
				guess = new_guess
				jacobian = jac(guess,na,nb,phi) #Evaluate the jacobian
				invjac = inv(jacobian) #Inverse jacobian
				f1 = fun(guess,na,nb,phi) #Calculate the function 
				new_guess = guess - damp*np.dot(invjac,f1) 

				#Tolerance criterion
				if abs(new_guess[0] - guess[0]) < 1e-8 and abs(new_guess[1]-guess[1]) < 1e-8: 
					x1[index] = phi
					x2[index] = new_guess[0]
					y2[index] = new_guess[1]
					break

		# Flips the function back
		#Numerical trick if na>nb
		if flipper ==1:
			x1 = 1 - x1
			x2 = 1 - x2

		n = np.size(x1) + 1
		x1 = np.reshape(np.append(x1,crit_phi),(n,1))
		x1=x1.tolist()
		x2=x2.tolist()

		#Has to reverse the order of x2, which was converted to a tuple in the previous line
		x2=x2[::-1] 

		#Adds crit chi to the end of y2
		y2 = np.reshape(np.append(y2,crit_chi),(n,1))
		y2=y2.tolist()
		y2i = y2[::-1]
		y2i.pop(0)

		#Concatenate the lists together
		phi = x1 + x2
		y2 = y2 + y2i


		return (phi,y2)

def fPlot(polya,polyb,na,nb,v0,jsondata):
	chi = float(request.form['chivalue'])
	chi = chi/v0
	type = jsondata['0']['type']

	fig = generate_figure(na,nb,chi)

	"""Add d3 stuff"""
	plugins.connect(fig, plugins.MousePosition())

	#Make the plot html/javascript friendly
	id = "fig01"
	json01 = json.dumps(mpld3.fig_to_dict(fig))

	list_of_plots = list()
	#Make a dictionary of plots
	plot_dict = dict()
	plot_dict['id'] = "fig01"
	plot_dict['json'] = json01
	list_of_plots.append(plot_dict)


	"""Spinodal"""
	crit_chi = .5*((1/(na**.5) + 1/(nb**.5))**2)
	if na == nb:
		crit_phi = 0.5
	else:
		crit_phi = (-nb + sqrt(na*nb) )/(na-nb)
	nav = 2./crit_chi
	nav = 1.0
	global flipper

	"""Flipper"""
	if na > nb:
			flipper = 1
			na, nb, w,x,y,z=  flip(na,nb,1,1,1,1)
	else:
			flipper = 0
	"""Set up the plot"""
	fig = plt.figure()
	fig.set_facecolor('white')
	axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
	axis.set_xlabel('Volume Fraction, \u03a6')
	axis.set_ylabel('Temperature, K')
	axis.set_title('Flory-Huggins Phase Diagram')

	"""Run Optimization"""
	x = np.arange(0.05,0.95,0.001)
	spinodal = nav*(.5*(1./(na*x) + 1./(nb-nb*x)))
	#spinodal = spinodal/nav

	if flipper == 1:
		x = 1 - x


	phi,y2 =  NR(na,nb,nav,crit_chi,flipper)
	#The line above and below do the same thing, one will replace the other soon.
	"""
	job = q.enqueue_call(
		func=NR, args=(na,nb,nav,crit_chi,flipper), result_ttl=5000)

	print job.get_id()
	"""

	#Convert list to np array
	y2 = np.asarray(y2)
	spinodal= np.asarray(spinodal)

	#Flip the plot w/ chi to be a function of temperature
	temp_unit = jsondata['0']['temperature_unit']

	if type == "Type 1":
		chi = float(jsondata['0']['chi'])

	
		y2 = (chi/v0)/y2
		spinodal = (chi/v0)/spinodal
		crit_chi = (chi/v0)/crit_chi

	elif type == "Type 2":
		y2 = float(jsondata['0']['chib']) / ( y2 - float(jsondata['0']['chia']))
		spinodal = float(jsondata['0']['chib']) / (spinodal - float(jsondata['0']['chia']))
		crit_chi= float(jsondata['0']['chib']) / (crit_chi- float(jsondata['0']['chia']))

	spinline = axis.plot(x,spinodal,'r',lw=2,label="Spinodal") 
	binline = axis.plot(phi,y2,'b',lw=2,label="Binodal")
	axis.legend()

	"""Add d3 stuff"""
	plugins.connect(fig, plugins.MousePosition())

	id2 = "fig02"
	json02 = json.dumps(mpld3.fig_to_dict(fig))

	#Attempt to make dictionary of plots
	plot_dict= dict()
	plot_dict['id'] = "fig02"
	plot_dict['json'] = json02
	list_of_plots.append(plot_dict)

	#Generate table
	zipped = zip(x,spinodal,y2)

	#Critical point
	critvals = [crit_phi,crit_chi]

	return jsondata,critvals,list_of_plots,zipped
