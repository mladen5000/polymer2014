from scipy.optimize import *
# Size of polymer
N = 1.
# chi parameter
w = 3.

# function defs
def freeenfun(x, w):
  phia = x[0]; phib = x[1]; phic = 1.-phia-phib
  return w*phia*phib + phia*log(phia) + phib/N*log(phib) + phic*log(phic)
def freeender(x, w):
  phia = x[0]; phib = x[1]; phic = 1.-phia-phib
  fda  = w*phib + log(phia/phic)
  fdb  = w*phia + log(phib)/N + 1./N - log(phic) - 1.
  return fda, fdb
def freeenderdir(t, x, y, w, lvl):
  phia = x[0] + t*cos(y); phib = x[1] + t*sin(y); phic = 1.-phia-phib
  fda, fdb = freeender(array([phia, phib]), w)
  fddir = fda*cos(y) + fdb*sin(y)
  return fddir - lvl
def freeenfun2(t, x, y, w):
  phia = x[0] + t*cos(y); phib = x[1] + t*sin(y); phic = 1.-phia-phib
  return freeenfun(array([phia, phib]), w)
def frompatophi(t, x, y):
  phia = x[0] + t*cos(y); phib = x[1] + t*sin(y); phic = 1.-phia-phib
  return phia, phib, phic

def getonepair(y, x, w): 
  y = y[0]
  a = cos(y)
  b = sin(y)
  bpol = roots(array([ (2*a**2*b**3 + 2*a**3*b**2)*w*N,
  ((4*a**2*b**2 + 2*a**3*b)*w*x[1] + (2*a*b**3 + 4*a**2*b**2)*w*x[0] - 2*a**2*b**2*w - a*b**3 - a**2*b**2)*N + a*b**3 + a**2*b**2,
  (2*a**2*b*w*x[1]**2 + ((4*a*b**2 + 4*a**2*b)*w*x[0] - 2*a**2*b*w - a*b**2)*x[1] + 2*a*b**2*w*x[0]**2 + (-2*a*b**2*w - b**3 - 2*a*b**2)*x[0] - a**2*b)*N + a*b**2*x[1] + (b**3 + 2*a*b**2)*x[0] - a*b**2,
  ((2*a*b*w*x[0] + a**2)*x[1]**2 + (2*a*b*w*x[0]**2 + (-2*a*b*w - b**2 - 2*a*b)*x[0] - a**2)*x[1])*N + b**2*x[0]*x[1] + b**2*x[0]**2 - b**2*x[0] ]))[1:]
  if any(iscomplex(bpol)):
    return array([nan, nan])
  lowerb = min(bpol)
  upperb = max(bpol)
  fdmax = freeenderdir(lowerb, x, y, w, 0)
  fdmin = freeenderdir(upperb, x, y, w, 0)
  hh = x[0]*tan(pi-y)
  ww = x[1]*tan(y-pi/2)
  if hh + x[1] < 1.:
    tmax = sqrt(x[0]**2 + hh**2)
  else:
    tmax = sin(pi/4.)/sin(3.*pi/4.-y)*(1-x[1]-x[0])
  if ww + x[0] < 1.:
    tmin = -sqrt(x[1]**2 + ww**2)
  else:
    tmin = -sin(pi/4.)/sin(y-3.*pi/4.)*(1-x[0]-x[1])
  def getonepair2(fd):
    ma = brentq(freeenderdir, tmin+1e-12, lowerb, args=(x, y, w, fd))
    mb = brentq(freeenderdir, upperb, tmax-1e-12, args=(x, y, w, fd))
    return array([ma, mb])
  def minfd(fd):
    onepair = getonepair2(fd)
    ypai0 = freeenfun2(onepair[0], x, y, w)
    ypai1 = freeenfun2(onepair[1], x, y, w)
    numd = (ypai1-ypai0)/diff(onepair)
    return (freeenderdir(onepair[0], x, y, w, 0) - numd)**2
  minres = minimize(minfd, (fdmax+fdmin)/2.)#, method='SLSQP', bounds=((fdmin, fdmax),))
  return getonepair2(minres.x)

def minang(x, w): 
  def minfun(y):
    onepair = getonepair(y, x, w)
    a = -onepair[0]/(onepair[1]-onepair[0])
    ff = (1-a)*freeenfun2(onepair[0], x, y, w) + a*freeenfun2(onepair[1], x, y, w)
    return ff
  minres = minimize(minfun, 135./180*pi)
  return minres.x, getonepair(minres.x, x, w)

aa = numpy.linspace(0.2, 0.4999, 256)
phia = ones((2, len(aa)))*0
phib = zeros_like(phia)
phic = zeros_like(phia)
for i in range(len(aa)):
  x = array([1, 1])*aa[i]; 
  y, tt = minang(x, w)
  phia[:, i], phib[:, i], phic[:, i] = frompatophi(tt, x, y)

kk = vstack((phia[0, :], phib[0, :], phic[0, :], phia[1, :], phib[1, :], phic[1, :])).transpose()
kk = kk[-isnan(kk[:, 0]), :]
savetxt('tern.dat', kk, header='A_a A_b A_c B_a B_b C_c', comments='')
