def EC_pdf_basic(x,y,x_obs,dx_obs,xtitle,ytitle) :
    import numpy as np
    from pylab import plot, legend, xlabel, ylabel, figure, savefig, title
    from pylab import mean, exp, sqrt, arange, ylim
    from lin_reg_DP import lin_reg_DP
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    import matplotlib
    import matplotlib.cm as cm
    from plot_scat_basic import plot_scat_basic
    import seaborn as sns
    from scipy.stats import norm
    
# Calculate mean and stdev of (equal model weight) prior
    mn_pr=mean(y)
    std_pr=np.std(y)
    print(mn_pr)
    print(std_pr)

    mn=x_obs
    std=dx_obs
    
# Calculate best-fit straight-line between x & y
    yf,a,b,da,db,xfit,yfit,yband=lin_reg_DP(x,y)
    
# Plot contours of probability fro linear regression
    mult=[-1,1]
    jconts=len(mult)
    figure()
    plot_scat_basic(x,y,xtitle,ytitle)
    plot(xfit,yfit,'k-.')
    plot(32,600,'b')
    for j in range (0,jconts):
      y1=yfit+mult[j]*yband
      plot(xfit,y1,'k--')
    pass

# Plot observational constraint
    xbest=x_obs
    xlo=x_obs-dx_obs
    xhi=x_obs+dx_obs
    plot([xbest,xbest],[300,900],'b-.')
    plot([xlo,xlo],[300,900],'b--',label='adJULES Constraint')
    plot([xhi,xhi],[300,900],'b--')
    legend(loc=3,prop={'size':11})
    savefig("Emergent_Relationship_Contours.pdf")

# Calculate PDF for IAV constraints
    x2=xfit
    nfitx=len(xfit)
    dx=x2[1]-x2[0]
    Px=x2
    Pi=3.142
    Px=1/sqrt(2*Pi*std**2) * exp(-((x2-mn)/(sqrt(2)*std))**2)

    miny=0.1*min(y)
    maxy=1.1*max(y)
    mfity=2000
    dy=(maxy-miny)/float(mfity)
    #y2=miny+dy*range(0,mfity)
    y2=miny+[dy* i for i in range(0,mfity)]

# Calculate "prior"
    Py_pr=y2
    Py_pr=1/sqrt(2*Pi*std_pr**2)*exp(-((y2-mn_pr)/(sqrt(2)*std_pr))**2)
    
# Calculate contours of probability in (x,y) space
    Pxy=np.zeros((nfitx,mfity))
    Pyx=np.zeros((mfity,nfitx))
    Py=np.zeros(mfity)
    Py_norm=0.0
    for m in range(0, mfity):
        Py[m]=0.0
        for n in range(0,nfitx):
            Py_given_x=1/sqrt(2*Pi*yband[n]**2) \
               * exp(-((y2[m]-yfit[n])/(sqrt(2)*yband[n]))**2)
            Pxy[n,m]=Px[n]*Py_given_x
            Pyx[m,n]=Pxy[n,m]
# Integrate over x to get Py
            Py[m]=Py[m]+Pxy[n,m]*dx
        pass
        Py_norm=Py_norm+Py[m]*dy
    pass

    print('NORMALISATION REQUIRED FOR POSTERIOR = ',Py_norm)
# Normalise Py
    for m in range(0, mfity):
        Py[m]=Py[m]/Py_norm
    pass

 #   figure()
 #   plot(x,y,'ro')
    z=100*Pyx
    CS = plt.contour(x2,y2,z)
    ylim(300,900) #ylim(min(y),max(y))
    plt.clabel(CS)
    savefig("Emergent_Relationship_contours_full.pdf")

# Plot PDF
    figure()
    plot(y2,Py,'k-',label='Emergent Constraint',linewidth=3)
    xlabel(ytitle,size=14)
    ylabel('Probablity Density',size=14)  
    dum=np.argmax(Py)
    ybest=y2[dum]
    a=np.column_stack((y2,Py))
    print(np.std(a))
    print(ybest-np.std(a))
    #plot(y2,norm.pdf(y2,ybest,91))
   
    dum_pr=np.argmax(Py_pr)
    ybest_pr=y2[dum_pr]
    #binny=min(y2)+(max(y2)-min(y2))*arange(11)/10.0
    binny=[150,250,350,450,550,650,750,850,950]
    binny=[100,200,300,400,500,600,700,800,900]
    n, bins, patches = plt.hist(y, bins=binny,\
                        normed=1, facecolor='orange',label='Equal-weight prior',edgecolor = "none")
    #sns.kdeplot(np.array(y),linewidth=3)
    legend(prop={'size':11})
    title("(a) PDF of Emergent Constraint",size=15,loc='left')
    savefig("Emergent_Constraint_PDF.pdf")  
 
# Calculate CDFs
    CDF=np.zeros(mfity)
    CDF_pr=np.zeros(mfity)
    CDF[0]=Py[0]*dy
    CDF_pr[0]=Py_pr[0]*dy
    for m in range(1,mfity):
        CDF[m]=CDF[m-1]+Py[m]*dy
        CDF_pr[m]=CDF_pr[m-1]+Py_pr[m]*dy
    pass

    
# Plot CDF   
    figure()
    plot(y2,CDF,'k-',label='Emergent Constraint',linewidth=3)
    
# Find 95% confidence limits    
    dum_up=CDF-0.975
    dum_975=dum_up**2
    dum_lo=CDF-0.025
    dum_025=dum_lo**2     
    val, n_lo = min((val, idx) for (idx, val) in enumerate(dum_025))
    val, n_hi = min((val, idx) for (idx, val) in enumerate(dum_975))
    val, n_best = max((val, idx) for (idx, val) in enumerate(Py))
    y_best=y2[n_best]
    y_lo=y2[n_lo]
    y_hi=y2[n_hi]
    
    plot([min(y2),max(y2)],[0.025,0.025],'k-.')
    plot([min(y2),max(y2)],[0.975,0.975],'k-.')
    xlabel(ytitle,size=14)
    ylabel('CDF (Probability<x)',size=14)
    n, bins, patches = plt.hist(y, bins=binny,cumulative=1,\
                                normed=1, facecolor='orange',\
                                label='Equal-weight prior',edgecolor = "none")
    legend(prop={'size':11})
    title("(b) CDF of Emergent Constraint",size=15,loc='left')
    savefig("Emergent_Constraint_CDF.pdf")
    
    print('EC on Y             = ',y_best,\
                          ',[',y_lo,'-',y_hi,']')


    return ybest, ybest_pr;

    
