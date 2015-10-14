#!/usr/bin/env python
import subprocess
import os

def hello_world():

	redis_url = os.getenv('REDISTOGO_URL','redis://localhost:6379')

	if redis_url == 'redis://localhost:6379':
		# This line is for local only
		"""
		output = subprocess.call(['./LOCAL_rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])
		
		"""
		var = 10 + 10
		output = str(var)
		return output
	
	else:
		print 'not local'
		print ''
		print ''

		#These 3 lines are for remote
		"""
		os.chdir('/app/SCFT_real')
		subprocess.call(['make'])
		output = subprocess.call(['./rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])
		"""
		output = 10 + 10
		output = str(output)


	return output
