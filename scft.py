#!/usr/bin/env python
import subprocess

def hello_world():

	subprocess.call(['make','-f','SCFT_real/Makefile'])
	output = subprocess.call(['./rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])
	#output = subprocess.check_output(['./rscft'])
	#print output
	#output = subprocess.check_output(['./SCFT_real/rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])

	""" Test Line """
	#output = subprocess.check_output('ls')
	return output
