import subprocess

def SCFT_execute():
	""" Make/Build the file onto heroku, and run the job"""
	print "FOUND IT"
	subprocess.call(['make','-f','SCFT_real/Makefile'])
	output = subprocess.call(['./rscft','37','3','6','3','3','outfile1','outfile2','infile','1.78'])
	return output
