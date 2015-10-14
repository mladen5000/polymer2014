#!/usr/bin/env python
import subprocess
import os
import time

import requests

def hello_world222():

	redis_url = os.getenv('REDISTOGO_URL','redis://localhost:6379')

	if redis_url == 'redis://localhost:6379':
		# This line is for local only
		"""
		output = subprocess.call(['./LOCAL_rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])
		
		"""
		time.sleep(12)
		var = 10 + 10
		output = str(var)
		return output
	
	else:
		#These 3 lines are for remote
		print 'notlocal'
		"""
		os.chdir('/app/SCFT_real')
		subprocess.call(['make'])
		output = subprocess.call(['./rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])
		"""
		output = 10 + 10


	return output
