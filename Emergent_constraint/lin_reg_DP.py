def lin_reg_DP(x,y) :
    import numpy as np
    from pylab import plot, legend, xlabel, ylabel
    from pylab import  figure, mean, exp, sqrt, sum
    from scipy import stats
    from lin_reg_DP import lin_reg_DP
    
# Based on "Least Squares fitting" equations from Math World website.
# This version checked against data on the Wikipedia "Simple Linear Regression" pages.
# It also calculates the +/- 1 sigma confidence limits in the regression [xfit,yband]
#
# IN THIS VERSION THE YBAND PREDICTION ERRORS ARE CALCULATED 
# ACCORDING TO DAVID PEARSON (AS USED IN COX ET AL., 2013)

    nx=len(x)
    ny=len(y)

    xm=mean(x)
    ym=mean(y)

    x2=x*x
    y2=y*y
    xy=x*y

    ssxx=sum(x2)-nx*xm*xm
    ssyy=sum(y2)-ny*ym*ym
    ssxy=sum(xy)-ny*xm*ym

    b=ssxy/ssxx
    a=ym-b*xm

    yf=a+b*x

    r2=ssxy**2/(ssxx*ssyy)
  
    e2=(y-yf)**2
    s2=sum(e2)/(nx-2)
  
    s=sqrt(s2)

    da=s*sqrt(1.0/nx+xm*xm/ssxx)
    db=s/sqrt(ssxx)


# Calculate confidence limits on fit (see Wikipedia page on "Simple Linear Regression")
#    minx=min(x)-1.1*(max(x)-min(x))
#    maxx=max(x)+1.1*(max(x)-min(x))
    minx=min(x)
    maxx=max(x)
    nfit=100
    dx=(maxx-minx)/nfit
    #xfit=minx+dx*range(0,nfit)
    xfit=minx+[dx*i for i in range(0,100)]
    #xfit=minx+dx*np.linspace(0,1,nfit-1)
    yfit=a+b*xfit
    yband=np.zeros(nfit)

# David Pearson's formula for "Prediction Error"
    for n in range (0,nfit):
        yband[n]=sqrt(s2*(1.0+1.0/nx+(xfit[n]-xm)**2/ssxx))
    pass

    return yf,a,b,da,db,xfit,yfit,yband;

