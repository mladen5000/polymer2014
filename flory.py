#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Numpy imports
from numpy.linalg import inv
from numpy import arange,asarray,zeros


#These can't be put __init__ for some reason
from SLCT import *
from VO import *
from FH import *
from structurefactor import structure_factor
import simpleA

#For talking to SQL
import json
import urllib

app = Flask(__name__)




@app.route('/')
@app.route('/index')
def index():
    return render_template("index.html")

@app.route('/saftdemo.html',methods=['POST','GET'])
def saftdemo():
	return render_template("saftdemo.html")

@app.route('/saftplot',methods=['POST','GET'])
def saftplot():
	if request.method == 'POST':

		#Initalization constants
		m = 2.457 #segment length
		sigma = 3.044 #segment diameter
		epsilon = 213.48 #well depth
		num_assocs = 3.0 #number of association sites
		ikappa = 1.0 #interaction strength paramater kappa
		ieps_ass = 1.0 #interaction strength paramater epsilon

		dens_num = 0.001

		#Put these into a single class
		compound = simpleA.Compound(sigma,epsilon,m,num_assocs,ikappa,ieps_ass,9999)



		#Set up demo
		eta = 0.05
		T = 0.2
		guess = [eta,T]

		#Generate Critical Point
		critvals = simpleA.findCrit(guess,dens_num,compound)
		Tc = critvals[1]
		Nc = critvals[0]

		#Generate Spinodal 
		Tvals, spin = simpleA.findSpin(Tc,Nc,dens_num,compound)
		Tvals =Tvals.tolist()
		spin =spin.tolist()

		while len(spin) > 100:
			del spin[::2]
			del Tvals[::2]


		#Set up figure and d3 plot 
		fig = Figure()
		fig.set_facecolor('white')
		axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
		axis.set_xlabel('Eta')
		axis.set_ylabel('Temperature, T')
		axis.set_title('SAFT LLEPhase Diagram')
		canvas = FigureCanvas(fig)

		#Plot spinodal
		spinline = axis.plot(spin,Tvals,'r',lw=2,label = "Spinodal")

		#Generate Binodal
		Tvals, bin = simpleA.findBinodal(dens_num,Tc,Nc,compound)
		Tvals =Tvals.tolist()
		bin =bin.tolist()

		while len(bin) > 200:
			del bin[::2]
			del Tvals[::2]

		#Plot Binodal
		binline = axis.plot(bin,Tvals,'b',lw=2, label = "Binodal") 

		axis.legend()
		#plugins.connect(fig, plugins.MousePosition())

		id1 = "fig01"
		json01 = json.dumps(mpld3.fig_to_dict(fig))

		#Attempt to make dictionary of plots
		plot_dict= dict()
		plot_dict['id'] = "fig01"
		pdid = "fig01"
		plot_dict['json'] = json01
		pdjson = json01
		

		#Generate table
		zipped2 = zip(Tvals,spin,bin)
			
		while len(zipped2) > 20:
			del zipped2[::2]

		#Critical point form
		critphi2 = critvals


		return render_template("exampleplots2.html",pdid=pdid,pdjson=pdjson,zipped2=zipped2,critphi2=critphi2)

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

		if request.form['eba'] == "":
			eba = 0
			ebb = 0

		else:
			eba = float(request.form['eba'])
			ebb = float(request.form['ebb'])


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

		#ADD
		#If checkmarked, take bending energies, else put in dummy value

		z = 6.0
		
		""" Parameters for specific polymers"""
		r1, p1,flex1 = SLCT_semiflex(polya,k1,m1,eba)
		r2, p2,flex2 =SLCT_semiflex(polyb,k2,m2,ebb)
		flex1.append(eba)
		flex2.append(ebb)

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
		if request.form['slctbutton'] == 'Generate Profile!':
			eps = float(request.form['eps'])

			fig = generate_SLCTfigure(na,nb,polya,polyb,k1,k2,m1,m2,eps)

			"""Add d3 stuff"""
			canvas = FigureCanvas(fig)
			output = StringIO.StringIO()
			canvas.print_png(output, bbox_inches='tight')
			plugins.connect(fig, plugins.MousePosition())


			id = "fig01"
			json01 = json.dumps(mpld3.fig_to_dict(fig))
			list_of_plots = list()

			#Attempt to make dictionary of plots
			plot_dict = dict()
			plot_dict['id'] = "fig01"
			plot_dict['json'] = json01
			list_of_plots.append(plot_dict)


			fig = Figure()
			fig.set_facecolor('white')
			axis = fig.add_subplot(1, 1, 1,axisbg='#f5f5f5')
			axis.set_xlabel('Volume Fraction, \u03a6')
			#axis.set_ylabel('Interaction Strength, \u03b5/kbT')
			axis.set_ylabel('Temperature, K')
			axis.set_title('SLCT Phase Diagram')

			print flex1, flex2
			"""Run Optimization"""
			phi, y2 = SLCT_NR(r1,r2,z,p1,p2,na,nb,flipper,eps,flex1,flex2)
		#	spinx,spiny = SLCT_Spinodal(r1,r2,z,p1,p2,na,nb,flipper)
			spinx,spiny = run_SLCT_flexspinodal(na,nb,flex1,flex2,eps,phi,flipper)
			

			"""Incorporate Epsilon"""
			#Convert list to np array
			y2 = np.asarray(y2)
			spiny= np.asarray(spiny)

			#Evaluate w/ epsilon
			y2 = eps/y2


			if flipper == 1:
				print "we gotta unflip!"
				na, nb, r1, r2, p1, p2 =  flip(na,nb,r1,r2,p1,p2)
				print r1, r2, z, p1, p2, na, nb

			"""Plot"""
			spinline = axis.plot(spinx,spiny,'r',lw=2,label="Spinodal") 
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
			plot_dict = dict()
			plot_dict['id'] = "fig02"
			plot_dict['json'] = json02
			list_of_plots.append(plot_dict)

			#Generate table
			zipped = zip(spinx,spiny,y2)

			#Critical point
			critphi = SLCT_crit(r1,r2,z,p1,p2,na,nb,eps)
			critphi = list(critphi)


			#return mpld3.fig_to_html(fig,template_type='simple')
			return render_template("slctplots.html",critphi=critphi,list_of_plots=list_of_plots,zipped=zipped)

@app.route('/plot', methods=['GET','POST'])	
def plot():
	#floryplot
	if request.method == 'POST':
		na = float(request.form['NFA'])
		nb = float(request.form['NFB'])

		#Talk to SQL
		#request the option selected from the dropdown
		polya = request.form['dropdowna']
		polyb = request.form['dropdownb']

		#create url string to fetch data from site
		webstr = "http://chidb.ci.uchicago.edu/basic_api.php?polymera="
		#had to switch because values are symmetrical
		site_url = webstr + polya + "&polymerb=" + polyb
		print site_url

		#open the site and load into json
		response = urllib.urlopen(site_url)
		jsondata = json.load(response)


		if request.form['florybutton'] == 'Generate Profile!':
			chi = float(request.form['chivalue'])
			type = jsondata['0']['type']

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
			nav = 1.0
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
			axis.set_ylabel('Temperature, K')
			axis.set_title('Flory-Huggins Phase Diagram')

			"""Run Optimization"""
			x = arange(0.05,0.95,0.001)
			spinodal = nav*(.5*(1./(na*x) + 1./(nb-nb*x)))
			#spinodal = spinodal/nav

			if flipper == 1:
				x = 1 - x


			phi,y2 =  NR(na,nb,nav,crit_chi,flipper)

			#Convert list to np array
			y2 = np.asarray(y2)
			spinodal= np.asarray(spinodal)

			#Flip the plot w/ chi to be a function of temperature
			temp_unit = jsondata['0']['temperature_unit']

			if type == "Type 1":
				chi = float(jsondata['0']['chi'])

				y2 = chi/y2
				spinodal = chi/spinodal
				crit_chi = chi/crit_chi

			elif type == "Type 2":
				y2 = float(jsondata['0']['chib']) / ( y2 - float(jsondata['0']['chia']))
				spinodal = float(jsondata['0']['chib']) / (spinodal - float(jsondata['0']['chia']))
				crit_chi= float(jsondata['0']['chib']) / (crit_chi- float(jsondata['0']['chia']))

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

			
			return render_template("exampleplots.html",polya=polya,polyb=polyb,jsondata=jsondata,critphi=critvals,list_of_plots=list_of_plots,zipped=zipped)

def flip(a1,a2,b1,b2,c1,c2):
		""" Switch a1 with a2, b1 with b2, c1 with c2"""
		return a2,a1,b2,b1,c2,c1

def vcrit():
	firstderiv = (1 + log(phi))/N + (-1)*(1 + log(1-phi-psi)) - (3*sigma*alpha/2.)*(sigma*phi + psi)**(1./2.)
	secondderiv = 1./(N*phi) + 1./(1-phi-psi) - (3*sigma*sigma*alpha/4.)*(sigma*phi + psi)**(-1./2.)
	thirdderiv = -1./(N*phi*phi) + (1./(1-phi-psi))**2 + (3*alpha*(sigma**3)/8.)*(sigma*phi + psi)**(-3./2.)

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
		h = np.zeros(( len(phi) ))
		s = np.zeros(( len(phi) ))
		g = np.zeros(( len(phi) ))

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
		r1, p1, r2, p2 = SLCT_constants(polya,polyb,k1,k2,m1,m2)
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

			return render_template("sfplot.html",id=id,json01=json01,zipped=zipped)


na = 1
nb = 1
crit_chi = .5 
crit_phi = 1
alpha = 3.655
N = 100



if __name__ == '__main__':
    app.run(debug=True)








