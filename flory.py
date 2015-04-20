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
from structurefactor import structure_factor


from flask import Flask,flash, request, make_response, render_template
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

app = Flask(__name__)
app.secret_key = 'mladen'



# Define a function that will return an HTML snippet.
def build_plot():
	#defined numpy.random as npr
    x_deets = npr.random(10)
    y_deets = npr.random(10)
    fig, ax = plt.subplots()
	#replaced pd with pandas, and imported pandas
    indata = pandas.DataFrame(x_deets, y_deets,)
    indata.plot(ax=ax)
	#replaced dumbs with json.dumps
    output = json.dumps(mpld3.fig_to_dict(fig))
    return output

# Define our URLs and pages.
@app.route('/blah')
def render_plot():
    sample_list = list(npr.randint(1,99999999,size=1))
    dict_of_plots = list()
    for i in sample_list:
        single_chart = dict()
        single_chart['id'] = "fig_" + str(i)
        print single_chart['id']

        single_chart['json'] = build_plot()     
        dict_of_plots.append(single_chart)
	id = "fig01"
	json01 = build_plot()
		#changed plots.html to exampleplots.html
    return render_template('exampleplots.html', id = id,json01=json01)#snippet=plot_snippet)

@app.route('/')
@app.route('/index')
def index():
    return render_template("index.html")

@app.route('/contact',methods=['GET','POST'])
def contact():
	form = ContactForm()
	if request.method == 'POST':
		if form.validate() == False:
			flash('FILL ALL FIELDS')
			return render_template('contact.html', form=form)
		else:
			return 'Form posted.'
	 
	elif request.method == 'GET':
		return render_template('contact.html', form=form)

	return render_template('contact.html', form=form)


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


		if request.form['vornbutton'] == 'Generate Profile!':
			psi = float(request.form['psi'])

			#I Use axis, instead of plt, same stuff though
			fig = Figure()
			fig.set_facecolor('white')
			axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
			axis.set_xlabel('Volume Fraction, \u03a6')
			axis.set_ylabel('Free Energy ')
			axis.set_title('Voorn-Overbeek  Diagram')
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

			id = "fig01"
			json01 = json.dumps(mpld3.fig_to_dict(fig))

			list_of_plots = list()
			#Attempt to make dictionary of plots
			plot_dict= dict()
			plot_dict['id'] = "fig01"
			plot_dict['json'] = json01
			list_of_plots.append(plot_dict)
			

			#return render_template("exampleplots.html",id=id,json01=json01)
			#return mpld3.show(fig)

		#elif request.form['vornbutton'] == 'Generate Phase!':
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
			
			#Critical point form
			critphi = vCriticalpoint(sigma,alpha,N)

			return render_template("exampleplots.html",critphi=critphi,list_of_plots=list_of_plots,zipped=zipped)

@app.route('/slct.html',methods=['POST','GET'])
def slct():
	return render_template("slct.html")

@app.route('/structurefactor.html',methods=['POST','GET'])
def structurefactor():
	return render_template("radiusgyration.html")


@app.route('/slctplot', methods=['GET','POST'])	
def slctplot():
	if request.method == 'POST':
		na = float(request.form['NFA'])
		nb = float(request.form['NFB'])
		polya = request.form['polya']
		polyb = request.form['polyb']
		print polya,polyb

		#Deal with k and m values for now
		if  polya =='PF' or polya == 'PG' or polya == 'PH' or polya == 'PI' or polya == 'PJ':
			k1 = float(request.form['k1'])
			m1 = float(request.form['m1'])
		else:
			k1 = 0
			m1 = 0

		if  polyb =='PF' or polyb == 'PG' or polyb == 'PH' or polyb == 'PI' or polyb == 'PJ':
			k2 = float(request.form['k2'])
			m2 = float(request.form['m2'])
		else:
			k2 = 0
			m2 = 0

		z = 6.0
		
		""" Parameters for specific polymers"""
		r1, p1, r2, p2 = SLCT_constants(polya,polyb,k1,k2)

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
			axis.set_ylabel('Interaction Strength, \u03b5/kbT')
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

		if request.form['florybutton'] == 'Generate Profile!':
			chi = float(request.form['chivalue'])

			fig = generate_figure(na,nb,chi)

			"""Add d3 stuff"""
			canvas = FigureCanvas(fig)
			output = StringIO.StringIO()
			canvas.print_png(output, bbox_inches='tight')
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
			axis.set_ylabel('Chi Parameter, \u03c7')
			axis.set_title('Flory-Huggins Phase Diagram')

			"""Run Optimization"""
			x = arange(0.05,0.95,0.001)
			spinodal = nav*(.5*(1./(na*x) + 1./(nb-nb*x)))
			spinodal = spinodal/nav

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
			
			return render_template("exampleplots.html",critphi=critvals,list_of_plots=list_of_plots,zipped=zipped)


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
		r1, p1, r2, p2 = SLCT_constants(polya,polyb,k1,k2)
		z = 6.0
		phi,h,s,g = SLCTfree(r1,r2,z,p1,p2,na,nb,eps)

		
		hline = axis.plot(phi,h,'r',lw=1,alpha = 0.5,label='Enthalpy')
		sline = axis.plot(phi,s,'b',lw=1,alpha = 0.5,label='Entropy')
		gline = axis.plot(phi,g,'g',lw=3,label="Free Energy")
		legend = axis.legend()
		return fig

@app.route('/sfplot', methods=['GET','POST'])	
def sfplot():
	if request.method == 'POST':
		na = float(request.form['Na'])
		nb = float(request.form['Nb'])
		ba = float(request.form['Ba'])
		bb = float(request.form['Bb'])
		phi = float(request.form['phia'])
		chi = float(request.form['chi'])


		if request.form['sfbutton'] == 'Generate Plot!':

			#I Use axis, instead of plt, same stuff though
			fig = Figure()
			fig.set_facecolor('white')
			axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
			axis.set_xlabel('q')
			axis.set_ylabel('S(q)')
			axis.set_title('Structure Factor')
			canvas = FigureCanvas(fig)
			q,sq = structure_factor(na,nb,ba,bb,phi,chi)

			
			gline = axis.plot(q,sq,'m',lw=3)
			legend = axis.legend()

			"""Add d3 stuff"""
			canvas = FigureCanvas(fig)
			output = StringIO.StringIO()
			canvas.print_png(output, bbox_inches='tight')
			plugins.connect(fig, plugins.MousePosition())

			id = "fig01"
			json01 = json.dumps(mpld3.fig_to_dict(fig))

			#Generate table
					
			zipped = zip(q,sq)

			while len(zipped) > 20:
				del zipped[::2]
				print zipped

			return render_template("sfplot.html",id=id,json01=json01,zipped=zipped)


na = 1
nb = 1
crit_chi = .5 
crit_phi = 1
alpha = 3.655
N = 100



if __name__ == '__main__':
    app.run(debug=True)








