

def vfun(x,phi1,sigma,alpha):
	"F1 = f'(phi_1a) - f'(phi_2a); F2 = (b-a)*f'(phi_1a) -[ f(phi_2a) - f(phi_1a) ]"
	return array([
		-phi1 + x[0] - m*(1 - phi1 - x[1]) + m*(1 - x[0] - x[1]) - (3.0/2)*alpha*m*sigma*sqrt(phi1*sigma + x[1]) + (1.0/2.0)*alpha*m*phi1*sigma*sqrt(phi1*sigma + x[1]) + 
  (1.0/2.0)*alpha*m*x[1]*sqrt(phi1*sigma + x[1]) + (3.0/2.0)*alpha*m*sigma*sqrt(sigma*x[0] + x[1]) - (1/2.0)*alpha*m*sigma*x[0]*sqrt(sigma*x[0] + x[1]) - 
  (1.0/2.0)*alpha*m*x[1]*sqrt(sigma*x[0] + x[1]) + log(phi1/2.0) - log(x[0]/2.0) + m*log(1 - phi1 - x[1]) - m*phi1*log(1 - phi1 - x[1]) - 
  m*(1 - phi1 - x[1])*log(1 - phi1 - x[1]) - m*x[1]*log(1 - phi1 - x[1]) - m*log(1 - x[0] - x[1]) + m*x[0]*log(1 - x[0] - x[1]) + 
  m*(1 - x[0] - x[1])*log[1 - x[0] - x[1]] + m*x[1]*log(1 - x[0] - x[1]),

	(phi1*log(phi1/2.0))/m - (log(x[0]/2.0)*x[0])/(1.0*m) + log[1 - phi1 - x[1]]*(1 - phi1 - x[1]) - log(1 - x[0] - x[1])*(1 - x[0] - x[1]) - 
  alpha*(phi1*sigma + x[1])**(3.0/2.0) + alpha*(sigma*x[0] + x[1])**(3.0/2.0) + 
  (-phi1 + x[0])*(-1 + 1.0/m + log(phi1/2.0)/(1.0*m) - log(1 - phi1 - x[1]) - (3.0/2.0)*alpha*sigma*sqrt(phi1*sigma + x[1]))
  ])

def vjac(x,alpha,N,phi1):
	"df1/dphi2, df1/dchi; df2/dphi2, df2/dchi"
	return array([[

	1.0 - 2*m - 1.0/x[0] + m*1.0/(1 - x[0] - x[1]) - (m*x[0])/(1 - x[0] - x[1]) - (m*x[1])/(1 - x[0] - x[1]) + 
  (3*alpha*m*sigma**2)/(4*sqrt(sigma*x[0] + x[1])) - (alpha*m*sigma**2*x[0])/(4*sqrt(sigma*x[0] + x[1])) - 
  (alpha*m*sigma*x[1])/(4*sqrt(sigma*x[0] + x[1])) - (1.0/2.0)*alpha*m*sigma*sqrt(sigma*x[0] + x[1]), # dF1/dphi2

	-(m*1.0/(1 - phi1 - x[1])) + (m*phi1)/(1 - phi1 - x[1]) + m*1.0/(1 - x[0] - x[1]) - (m*x[0])/(1 - x[0] - x[1]) + (m*x[1])/(1 - phi1 - x[1]) - 
  (m*x[1])/(1 - x[0] - x[1]) - (3*alpha*m*sigma)/(4*sqrt(phi1*sigma + x[1])) + (alpha*m*phi1*sigma)/(4*sqrt(phi1*sigma + x[1])) + 
  (alpha*m*x[1])/(4*sqrt(phi1*sigma + x[1])) + (1.0/2.0)*alpha*m*sqrt(phi1*sigma + x[1]) + (3*alpha*m*sigma)/(4*sqrt(sigma*x[0] + x[1])) - 
  (alpha*m*sigma*x[0])/(4*sqrt(sigma*x[0] + x[1])) - (alpha*m*x[1])/(4*sqrt(sigma*x[0] + x[1])) - (1.0/2.0)*alpha*m*sqrt(sigma*x[0] + x[1])
		 	], #dF1/dpsi

	[log(phi1/2.0)/m - log(x[0]/2.0)/(m*1.0) - log(1 - phi1 - x[1]) + log(1 - x[0] - x[1]) - (3.0/2.0)*alpha*sigma*sqrt(phi1*sigma + x[1]) + 
  (3.0/2.0)*alpha*sigma*sqrt(sigma*x[0] + x[1]) , #dF2/dphi2

	log(phi1/2.0)/m - log(x[0]/2.0)/m - log(1 - phi1 - x[1]) + log(1 - x[0] - x[1]) - (3.0/2.0)*alpha*sigma*sqrt(phi1*sigma + x[1]) + 
  (3.0/2.0)*alpha*sigma*sqrt(sigma*x[0] + x[1])

			]]) #dF2/dpsi

