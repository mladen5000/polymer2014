#!/usr/bin/env python
# -*- coding: utf-8 -*-


#Numpy imports
from numpy.linalg import inv
from numpy import arange,asarray,zeros


#These can't be put __init__ for some reason
import SLCT
import VO
from FH import *
from structurefactor import structure_factor
import simpleA

#Task queueing, redis
from rq import Queue
from rq.job import Job
from worker import conn
from rq_dashboard import RQDashboard

#For talking to chiSQL
import json
import urllib

#Used for Self Consistent Field Theory
import subprocess
from scft import *

#Used to Debug
from flask_debugtoolbar import DebugToolbarExtension
import saftdemocode

#Extra
from general_route_functions import *



#Initialize
app = Flask(__name__)
#redis queue, conn is the name of the redis connection as defined on worker
q = Queue(connection=conn)
app.config['REDIS_URL'] = os.getenv('REDISTOGO_URL', 'redis://localhost:6379')
app.config['RQ_POLL_INTERVAL'] = 5
app.config['DEBUG'] = True
app.config['SECRET_KEY'] = 'use a better key'

RQDashboard(app,'/rq')
toolbar = DebugToolbarExtension(app)





@app.route('/')
@app.route('/index')
def index():
    return render_template("index.html")

@app.route('/scft_enqueue',methods=['GET','POST'])
def scft_enqueue():
	if request.method == 'POST':
		if request.form['enqueuebutton'] == 'Submit':
			print 'queueing!'
			job = q.enqueue_call(
				func=hello_world, args=(), result_ttl=5000,timeout=9999
			)
			return render_template("scft_enqueue.html")
		else:
			print 0

	else:
		print 'NOT POST'
		return render_template("scft_enqueue.html")


@app.route('/rs/<job_key>',methods=['GET'])
def get_results(job_key):
	print job_key
	job = Job.fetch(job_key,connection=conn)

	if job.is_finished:
		return job.result
	else:
		return "Nay!"


### SAFT DEMO ###
@app.route('/saftdemo.html',methods=['POST','GET'])
def saftdemo():
	return render_template("saftdemo.html")

@app.route('/saftplot',methods=['POST','GET'])
def saftplot():
	if request.method == 'POST':
		#Call the function in saftdemocode.py
		pdid,pdjson,zipped2,critphi2 = saftdemocode.saft_plot()

		return render_template("exampleplots2.html",pdid=pdid,pdjson=pdjson,zipped2=zipped2,critphi2=critphi2)

#################

### GENERAL HTML PAGES ###
@app.route('/howto.html')
def howto():
	return render_template("howto.html")

@app.route('/flory.html',methods=['POST','GET'])
def flory():
	return render_template("flory.html")

@app.route('/apps.html')
def apps():
	return render_template("apps.html")

#####################


### VORN PLOT ####
@app.route('/vorn.html',methods=['POST','GET'])
def vorn():
	return render_template("vorn.html")

@app.route('/vornplot', methods=['GET','POST'])	
def vornplot():
	if request.method == 'POST':
		N = float(request.form['N'])
		sigma = float(request.form['sigma'])
		psi = float(request.form['psi'])
		if request.form['vornbutton'] == 'Generate Profile!':

			critphi, list_of_plots, zipped = VO.vPlot(3.655,sigma,psi,N)
			return render_template("slctplots.html",critphi=critphi,list_of_plots=list_of_plots,zipped=zipped)

###################


### SIMPLE LATTICE CLUSTER, SLCT ###
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
		r1, p1,flex1 = SLCT.SLCT_semiflex(polya,k1,m1,eba)
		r2, p2,flex2 =SLCT.SLCT_semiflex(polyb,k2,m2,ebb)
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

#######################################


@app.route('/plot', methods=['GET','POST'])	
def plot():
	#floryplot
	if request.method == 'POST':
		na = float(request.form['NFA'])
		nb = float(request.form['NFB'])
		v0 = float(request.form['v0'])

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
			chi = chi/v0
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
			#The line above and below do the same thing, one will replace the other soon.
			job = q.enqueue_call(
				func=NR, args=(na,nb,nav,crit_chi,flipper), result_ttl=5000)

			print job.get_id()

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


@app.route("/results/<job_key>", methods=['GET'])
def fget_results(job_key):
	#workerfloryplot
    job = Job.fetch(job_key, connection=conn)

    if job.is_finished:
        return jsonify(job.result), 200
    else:
        return "Nay!", 202

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
    app.run()








