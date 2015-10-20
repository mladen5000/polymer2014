#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Flask
from flask import Flask, request, render_template

#Numpy imports
from numpy.linalg import inv
from numpy import arange,asarray,zeros

#Modules
import models.SLCT 
import models.VO as VO
from models.FH import *
from models.structurefactor import structure_factor
import models.simpleA
from models.scft import *
from models.general_route_functions import *
import models.saftdemocode

#Task queueing, redis
from rq import Queue
from rq.job import Job
from worker import conn
from rq_dashboard import RQDashboard
import subprocess

#For talking to chiSQL
import json
import urllib

#Used to Debug
from flask_debugtoolbar import DebugToolbarExtension



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

		critphi, list_of_plots, zipped =SLCT.sPlot(polya,
					polyb,na,nb)

		return render_template("slctplots.html",
				critphi=critphi,list_of_plots=list_of_plots,zipped=zipped)

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

			jsondata,critvals,list_of_plots,zipped = fPlot(polya,polyb,na,nb,v0,jsondata)
			
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








