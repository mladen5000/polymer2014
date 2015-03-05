#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import random
from math import *
from numpy import *
from numpy.linalg import inv
import matplotlib.pyplot as plt
import StringIO
import mpld3
from mpld3 import plugins

from SLCT import *
from VO import *
from FH import *


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

@app.route('/vornplot', methods=['GET','POST'])	
def vornplot():
	if request.method == 'POST':
		N = float(request.form['N'])
		sigma = float(request.form['sigma'])


		if request.form['vornbutton'] == 'Generate Free Energy!':
			
			psi = float(request.form['psi'])

			#I Use axis, instead of plt, same stuff though
			fig = Figure()
			fig.set_facecolor('white')
			axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
			axis.set_xlabel('Volume Fraction, \u03a6')
			axis.set_ylabel('Free Energy ')
			axis.set_title('Voorn-Overbeek Phase Diagram')
			canvas = FigureCanvas(fig)
			phi,h,s,g = vfree(N,psi,sigma)

			
			hline = axis.plot(phi,h,'r',lw=1,alpha = 0.5,label='Enthalpy')
			sline = axis.plot(phi,s,'b',lw=1,alpha = 0.5,label='Entropy')
			gline = axis.plot(phi,g,'g',lw=3,label="Free Energy")
			legend = axis.legend()

			"""Add d3 stuff"""
			canvas = FigureCanvas(fig)
			output = StringIO.StringIO()
			canvas.print_png(output, bbox_inches='tight')
			plugins.connect(fig, plugins.MousePosition())

			return mpld3.fig_to_html(fig,template_type='simple')

		elif request.form['vornbutton'] == 'Generate Phase!':
			"""Set up the plot"""
			fig = Figure()
			fig.set_facecolor('white')
			axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
			axis.set_xlabel('Volume Fraction, \u03a6')
			axis.set_ylabel('Salt Concentration, \u03a8')
			axis.set_title('Voorn-Overbeek Phase Diagram')
			canvas = FigureCanvas(fig)

			"""Move Spinodal Elsewhere"""
			phi,y2 =  vNR(alpha,N,sigma)
			x, spinodal = vSpinodal(sigma,alpha,N)
			spinline = axis.plot(x,spinodal,'r',lw=2,label = "Spinodal")
			binline = axis.plot(phi,y2,'b',lw=2, label = "Binodal") 
			axis.legend()

			"""Make this organized like the other stuff"""
			plugins.connect(fig, plugins.MousePosition())
			return mpld3.fig_to_html(fig)

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

		z = 6.0
		
		""" Parameters for specific polymers"""
		r1, p1, r2, p2 = SLCT_constants(polya,polyb)

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
		if request.form['slctbutton'] == 'Generate Free Energy!':
			eps = float(request.form['eps'])

			fig = generate_SLCTfigure(na,nb,polya,polyb,k1,k2,m1,m2,eps)

			"""Add d3 stuff"""
			canvas = FigureCanvas(fig)
			output = StringIO.StringIO()
			canvas.print_png(output, bbox_inches='tight')
			plugins.connect(fig, plugins.MousePosition())
			return mpld3.fig_to_html(fig,template_type='simple')

		else:
			"ELSE PART HAPPENED"

		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Volume Fraction, \u03a6')
		axis.set_ylabel('\u03b5/kbT')
		axis.set_title('SLCT Phase Diagram')

		"""Run Optimization"""
		phi, y2 = SLCT_NR(r1,r2,z,p1,p2,na,nb,flipper)
		spinx,spiny = SLCT_Spinodal(r1,r2,z,p1,p2,na,nb,flipper)


		"""Plot"""
		spinline = axis.plot(spinx,spiny,'r',lw=2,label="Spinodal") 
		binline = axis.plot(phi,y2,'b',lw=2,label="Binodal")
		axis.legend()

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

		if request.form['florybutton'] == 'Generate Free Energy!':
			chi = float(request.form['chivalue'])

			fig = generate_figure(na,nb,chi)

			"""Add d3 stuff"""
			canvas = FigureCanvas(fig)
			output = StringIO.StringIO()
			canvas.print_png(output, bbox_inches='tight')
			plugins.connect(fig, plugins.MousePosition())

			return mpld3.fig_to_html(fig,template_type='simple')

		elif request.form['florybutton'] == 'Generate Phase!':

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
			axis.set_xlabel('Volume Fraction, \u03a6')
			axis.set_ylabel('\u03c7')
			axis.set_title('Flory-Huggins Phase Diagram')

			"""Run Optimization"""
			x = arange(0.05,0.95,0.001)
			spinodal = nav*(.5*(1./(na*x) + 1./(nb-nb*x)))

			if flipper == 1:
				x = 1 - x

			phi,y2 =  NR(na,nb,nav,crit_chi,flipper)
			spinline = axis.plot(x,spinodal,'r',lw=2,label="Spinodal") 
			binline = axis.plot(phi,y2,'b',lw=2,label="Binodal")
			axis.legend()

			"""Add d3 stuff"""
			canvas = FigureCanvas(fig)
			output = StringIO.StringIO()
			canvas.print_png(output, bbox_inches='tight')
			plugins.connect(fig, plugins.MousePosition())

			return mpld3.fig_to_html(fig,template_type='simple')


def flip(a1,a2,b1,b2,c1,c2):
		""" Switch a1 with a2, b1 with b2, c1 with c2"""
		return a2,a1,b2,b1,c2,c1




def vcrit():
	firstderiv = (1 + log(phi))/N + (-1)*(1 + log(1-phi-psi)) - (3*sigma*alpha/2.)*(sigma*phi + psi)**(1./2.)
	secondderiv = 1./(N*phi) + 1./(1-phi-psi) - (3*sigma*sigma*alpha/4.)*(sigma*phi + psi)**(-1./2.)
	thirdderiv = -1./(N*phi*phi) + (1./(1-phi-psi))**2 + (3*alpha*(sigma**3)/8.)*(sigma*phi + psi)**(-3./2.)

"Flory huggins formula"
def flory_G(phi,na,nb,chi):
	"""Plots free energy"""
	enthalpy = chi*phi*(1-phi)
	entropy = phi/na * log(phi) + (1.-phi)/nb * log(1-phi) 
	f = phi/na * log(phi) + (1.-phi)/nb * log(1-phi) + chi*phi*(1-phi)
	return enthalpy,entropy,f

def generate_figure(na,nb,chi):
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Volume Fraction, \u03a6')
		axis.set_ylabel('Free Energy')
		axis.set_title('Flory-Huggins Free Energy Diagram')

		"""Run Optimization"""
		"""Need to move these lines it's own function"""
		phi = arange(0.0001,0.99,0.001)
		h = zeros(( len(phi) ))
		s = zeros(( len(phi) ))
		g = zeros(( len(phi) ))

		i = 0
		for current_phi in phi:
			h[i],s[i],g[i] = flory_G(current_phi,na,nb,chi)
			i += 1


		hline = axis.plot(phi,h,'r',lw=1,alpha = 0.5,label='Enthalpy')
		sline = axis.plot(phi,s,'b',lw=1,alpha = 0.5,label='Entropy')
		gline = axis.plot(phi,g,'g',lw=3,label="Free Energy")
		legend = axis.legend()
		return fig

def generate_SLCTfigure(NFA,NFB,polya,polyb,k1,k2,m1,m2,eps):
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Volume Fraction, \u03a6')
		axis.set_ylabel('Free Energy')
		axis.set_title('SLCT Free Energy Diagram')

		"""Run Optimization"""
		"""Need to move these lines it's own function"""
		r1, p1, r2, p2 = SLCT_constants(polya,polyb)
		z = 6.0
		phi,h,s,g = SLCTfree(r1,r2,z,p1,p2,na,nb,eps)

		
		hline = axis.plot(phi,h,'r',lw=1,alpha = 0.5,label='Enthalpy')
		sline = axis.plot(phi,s,'b',lw=1,alpha = 0.5,label='Entropy')
		gline = axis.plot(phi,g,'g',lw=3,label="Free Energy")
		legend = axis.legend()
		return fig

na = 1
nb = 1
crit_chi = .5 
crit_phi = 1
alpha = 3.655
N = 100



if __name__ == '__main__':
    app.run(debug=True)








