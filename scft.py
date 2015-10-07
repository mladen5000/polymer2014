#!/usr/bin/env python
import subprocess

def hello_world():
	#open a file to write, execute the C code and stdout it
	#subprocess.call(['make','-f','SCFT_real/local_make'])
	"""
	with open('out-file.txt', 'w') as f:
	    subprocess.call(['./SCFT_real/rscft','37','1.78','1.78','1.78','0','outfile1','outfile2','infile','1.78'], stdout=f)
	f = open('out-file.txt','r')
	"""
	#output = subprocess.check_output(['./SCFT_real/rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])
	#subprocess.call(['rm','SCFT_real/rscft','SCFT_real/rscft.o', 'SCFT_real/main.o'])
	output = subprocess.check_output('ls')
	return output
