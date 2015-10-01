#!/usr/bin/env python
from math import *

def function1(n,Ebb,Ebs):
	#u2b is defined as N2b/M, the number of 2-bonds on the backbone
	#u2s is defined as N2s/M, the number of 2-bonds on the sidechain 
	#u2d is defined as N2d/M, the number of 2-bonds where 1 is on the sidechain and 1 is on the backbone
	#M is the total number of bonds in a chain

	#Going to use constants for now, simple example for poly(n-alpha olefins)

	#Energy, Eb, and epsilon, eps, are given in terms of Eb/kb and eps/kb
	Ebb = Ebb
	Ebs = Ebs

	#n, number of UAG (lattice sites occupied by side chain
	n = 4.
	z = 6.0

	#p, polymerization index
	p = 1.0

	T = 300.0


	#Counting indices, these are valid for poly(n-alpha olefins)
	M = (n+2)*p + 1
	N1 = M - 1
	print M
	N2s = (n-1)*p
	N2b = 2*p - 1
	N2d = 2*p

	#These get used
	Zbb = 1.0 + (z/2.0 - 1.0)*exp(-Ebb/(T)) 
	Zbs = 1.0 + (z/2.0 - 1.0)*exp(-Ebs/(T)) 

	print "ZBB,ZBS"
	print Zbb,Zbs

	#Calculate gi
	gb = (z/2)*exp(Ebb/T)/Zbb
	gs = (z/2)*exp(Ebs/T)/Zbs
	print "GB,GS"
	print gb,gs


	#These get used
	u1 = N1/M
	u2b = N2b/M
	u2s = N2s/M
	u2d = N2d/M

	u2t = gb*u2b + gs*u2s + u2d
	print "u1",u1
	print "u2b",u2b
	print "u2s",u2s
	print "u2d",u2d
	print "u2t",u2t
	


	#Calculate F0
	sumlogs = u2b*log(Zbb) + u2s*log(Zbs)
	f0_1 = (phi/M)*log( 2*phi / (M*z**(p+1)) ) 
	f0_2 = phi - (phi/M) + (1-phi)*log(1-phi)
	f0_3 = (z/2)*beta*eps*phi**2 - sumlogs
	f0 = f0_1 + f0_2 + f0_3
	
	#Calculate F1
	f1_1 = (-1*u2t*phi + (u1*phi)**2)/z
	f1_2 = u1*beta*eps*phi*(1-phi)**2
	f1_3 = z*(1-phi)**2*(beta*eps*phi)**2 / 4.0
	f1 = f1_1 + f1_2 + f1_3

	#Free energy, F / (NlkT)

	free_energy = f0 - f1
	
function1(2,400.0,800.0)


