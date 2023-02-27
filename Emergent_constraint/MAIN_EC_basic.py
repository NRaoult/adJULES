#!/usr/bin/env python
import numpy as np
from pylab import plot, show, bar, legend, colors, axes, xlabel, ylabel, hist
from pylab import title, savefig, axis, figure, semilogx, mean, exp, sqrt
from pylab import log, arctan
from scipy import stats
from plot_scat_basic import plot_scat_basic
from EC_pdf_basic import EC_pdf_basic

# Read x and y variables
fname='Booth_2012.csv'
xy_data=np.loadtxt(fname, \
                 dtype={'names': ('xcol','ycol'),'formats' : ('f8','f8')},\
                 comments='#',delimiter=',', converters=None, skiprows=1, \
                 usecols=None, unpack=False, ndmin=0)
xlab='T$_\mathregular{opt}$ (K)'
ylab='CO$_2$ Change, 2100-1900 (ppmv)'
x=xy_data['xcol']
y=xy_data['ycol']
npts=len(x)

# Input obsrvational estimate of x
x_obs=35.392 #35
dx_obs=0.922

#x_obs=43.027
#dx_obs=3.376

#x_obs=32
#dx_obs = sqrt(25/sqrt(0.5))

#x_obs=33.649
#dx_obs=0.628

#x_obs=37.17861#
#dx_obs= 2.897104

#x_obs=37.35398
#dx_obs=1.24051
# Plot Emergent Relationship
figure()
plot_scat_basic(x,y,xlab,ylab)
deg=1
p=np.polyfit(x,y,deg)
yfit=np.zeros(npts)
for n in range(0,npts):
    yfit[n]=0.0
    for j in range(0,deg+1):
        yfit[n]=yfit[n]+p[deg-j]*x[n]**j
    pass
pass
plot(x,yfit,'k-')
xbest=x_obs
xlo=x_obs-dx_obs
xhi=x_obs+dx_obs
plot([xbest,xbest],[min(y),max(y)],'b-.')
plot([xlo,xlo],[min(y),max(y)],'b--',label='ADJULES Constraint')
plot([xhi,xhi],[min(y),max(y)],'b--')
legend(loc=3,prop={'size':11})
savefig("Emergent_Relationship.pdf")

# Calculate PDF from Emergent Constraint
ytit=ylab
xtit=xlab
ybest,ybest_pr=EC_pdf_basic(x,y,x_obs,dx_obs,xlab,ylab)
print(' Equal-weight prior Estimate of Y = ',ybest_pr)
print('Emergent Constraint Estimate of Y = ',ybest)
