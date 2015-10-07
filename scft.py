#!/usr/bin/env python
import subprocess
import os

def hello_world():


	#These 3 lines are for remote
	os.chdir('/app/SCFT_real')
	subprocess.call(['make'])
	output = subprocess.call(['./rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])

	# This line is for local only
	#output = subprocess.call(['./LOCAL_rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])

	return output
