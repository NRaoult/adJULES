def plot_scat_basic(x,y,xlab,ylab) :
    import numpy as np
    from pylab import plot, show, bar, legend, colors, axes, xlabel, ylabel
    from pylab import title, savefig, axis, figure, semilogx, mean, exp, sqrt
    from pylab import log, arctan
    import scipy
    import matplotlib.pyplot as plt
         
   
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot(x,y,'ro')
    xlabel(xlab,size=14)
    ylabel(ylab,size=14)
       
    return;

