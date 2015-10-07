#!/usr/bin/env python
import subprocess
import os

def hello_world():


	#os.chdir('/app/SCFT_real')
	#subprocess.call(['make'])
	output = subprocess.call(['./LOCAL_rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])
	#output = subprocess.check_output(['./rscft'])
	#print output
	#output = subprocess.check_output(['./SCFT_real/rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])

	""" Test Line """
	#output = subprocess.check_output('ls')
	return output
