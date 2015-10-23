#!/usr/bin/env python

from math import *
from numpy import *
from numpy.linalg import inv
import matplotlib.pyplot as plt
import StringIO
import mpld3
from mpld3 import plugins

from flask import Flask, Response, request, redirect, make_response, render_template, jsonify
from sqlalchemy import create_engine

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import json
import os
#import MySQLdb

#app = Flask(__name__)
app = Flask(__name__, static_url_path='')

# Database is currently local to my machine, won't work for everyone
DEBUG = 1


if DEBUG == 0:
	db = MySQLdb.connect(host="localhost", user="roselyne", db="chidb", passwd="roselyne", use_unicode=True, charset="utf8")
	cur = db.cursor()

# Root path
@app.route('/')
def root():
	return render_template('index.html')

# Polymer page
@app.route('/polymer', methods=['POST','GET'])
def polymer():
	if request.method == 'POST':
     		poly = request.form.get('poly')
     	htmlpath = "/polymer.html?poly="+poly
     	return redirect(htmlpath)

# Get chi values for polymer
@app.route('/chis',methods=['POST','GET'])
def chis():
	if request.method == 'POST':
		if DEBUG == 0:
			poly = request.get_json()
			poly = poly['polymer']
			query0 = "Select compound1, compound2, chinumber, chia, chib, chic, method, notes, doi from reviewed_chis where (compound1 = '%s' or compound2='%s' )and type1='polymer' and type2='polymer'" % (poly, poly)
			cur.execute(query0)
			db.commit()
		        results0 = cur.fetchall()

			# Form table
		    	tablecontent = '<table id="chitable" class="display" cellspacing="0" width="100%">'
	    		chi = u'\u03C7' 
	    		tableheader = "<thead><tr><th>Polymer</th><th>%s</th><th>%sa</th><th>%sb</th><th>%sc</th><th>Method</th><th>Note</th><th>DOI</th></tr></thead>" % (chi, chi, chi, chi)
	    
	    		tablecontent = tablecontent + tableheader
	    		for row0 in results0:
				# Get link to paper from second table
				query1 = "select remotepath from papers where doi='%s'" % (row0[8])
				cur.execute(query1)
				db.commit()
				results1 = cur.fetchone()	
				url = results1[0]

				# Form row
				tablerow = "<tr>"
	
				# Get other half of the polymer pair
				if row0[0].lower().find(poly) == -1:
					polymer = row0[0]
				else:
					polymer = row0[1]
				tablerow = tablerow + "<td>" + polymer + "</td>"
				tablerow = tablerow + "<td>" + str(row0[2]) + "</td>"
				tablerow = tablerow + "<td>" + str(row0[3]) + "</td>"
				tablerow = tablerow + "<td>" + str(row0[4]) + "</td>"
				tablerow = tablerow + "<td>" + str(row0[5]) + "</td>"
				tablerow = tablerow + "<td>" + str(row0[6]) + "</td>"
				tablerow = tablerow + "<td>" + row0[7] + "</td>"
				tablerow = tablerow + "<td><a target='_blank' href='"+url+"'>" + row0[8] + "</a></td>"
				tablerow = tablerow + "</tr>"
				tablecontent = tablecontent + tablerow	
	    		tablecontent = tablecontent + "</table>"
	   	else: 
			tablecontent = "<h4>DEBUG on - DB not connected: Showing results for butadiene </h4></br><table id=\"chitable\" class=\"display\" cellspacing=\"0\" width=\"100%\"><thead><tr><th>Polymer</th><th>X</th><th>Xa</th><th>Xb</th><th>Xc</th><th>Method</th><th>Note</th><th>DOI</th></tr></thead><tr><td>polystyrene</td><td>0.047</td><td>0.0</td><td>0.0</td><td>0.0</td><td></td><td></td><td><a target='_blank' href='http://pubs.acs.org/doi/full/10.1021/ma401957s'>10.1021/ma401957s</a></td></tr><tr><td>styrene</td><td>25.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>scattering1</td><td>X is given in terms of XN</td><td><a target='_blank' href='http://pubs.acs.org/doi/full/10.1021/ma3004136'>10.1021/ma3004136</a></td></tr><tr><td>polystyrene</td><td>0.05</td><td>0.0</td><td>0.0</td><td>0.0</td><td>scattering1</td><td>at 0 mol % C6F13H modifications</td><td><a target='_blank' href='http://pubs.acs.org/doi/full/10.1021/ma400533w'>10.1021/ma400533w</a></td></tr><tr><td>polystyrene</td><td>0.75</td><td>0.0</td><td>0.0</td><td>0.0</td><td>scattering1</td><td>at 80 mol % C6F13H modification</td><td><a target='_blank' href='http://pubs.acs.org/doi/full/10.1021/ma400533w'>10.1021/ma400533w</a></td></tr><tr><td>polystyrene</td><td>0.0555</td><td>0.0</td><td>0.0</td><td>0.0</td><td>scattering1</td><td>Composition given in volume fractions. Temperature cited as room temperature, which is 23 degrees Celsius on average.</td><td><a target='_blank' href='http://pubs.acs.org/doi/full/10.1021/ma500633b'>10.1021/ma500633b</a></td></tr><tr><td>polystyrene</td><td>0.16</td><td>0.0</td><td>0.0</td><td>0.0</td><td>md_simulation</td><td></td><td><a target='_blank' href='http://pubs.acs.org/doi/full/10.1021/ma202717p'>10.1021/ma202717p</a></td></tr><tr><td>polystyrene</td><td>32.5</td><td>0.0</td><td>0.0</td><td>0.0</td><td>micro_phase_tri</td><td></td><td><a target='_blank' href='http://pubs.acs.org/doi/full/10.1021/ma301487e'>10.1021/ma301487e</a></td></tr></table>"
	return tablecontent	

	


if __name__ == '__main__':
    app.run(debug=True)
















