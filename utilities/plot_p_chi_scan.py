import numpy as np
from handle_netcdf_p_chi_grid import load_p_chi_grid
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
plt.rcParams['savefig.dpi'] = 75

plt.rcParams['figure.autolayout'] = False
plt.rcParams['figure.figsize'] = 10, 6
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['lines.markersize'] = 8
plt.rcParams['legend.fontsize'] = 14

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"
plt.rcParams['font.serif'] = "cm"
plt.rcParams['text.latex.preamble'] = "\usepackage{subdepth}, \usepackage{type1cm}"


def plotFigs_chi_to_p(fnames,labels,nostochdata,save=None):
    nostochdata=np.loadtxt(nostochdata).T
    plt.figure(figsize=(5*1.62,5))
    plt.fill_betweenx(nostochdata[0],nostochdata[1],np.amax(nostochdata[1]),color='g',label=r'$1:1$ response')
    plt.fill_betweenx(nostochdata[0],nostochdata[2],nostochdata[1],color='c',label=r'$2:1$ response')
    plt.fill_betweenx(nostochdata[0],nostochdata[3],nostochdata[2],color='m',label=r'Chaotic response')
    plt.plot(nostochdata[3],nostochdata[0],'b',label=r'Periodic collapse')
    cl = plt.cm.autumn_r(np.linspace(0, 1, len(fnames)))
    for i,fname in enumerate(fnames):
        if fname.endswith(".nc"):
            fname=fname[:-3]
        data,p,chi=load_p_chi_grid(fname+".nc")
        p_stoch_col=p[np.argmax(data,axis=0)]
        p_stoch_col_new = interp1d(chi,p_stoch_col, kind='quadratic')(nostochdata[0])
        plt.plot(p_stoch_col_new,nostochdata[0],'r',color=cl[i],label=r'$\sigma={:3.2f}$'.format(labels[i]))
    plt.ylim([0,1])
    plt.xlim([np.amin(nostochdata[3]),np.amax(nostochdata[1])])
    plt.xlabel(r'$p$')
    plt.ylabel(r'$\chi$')
    plt.legend(loc='lower left')
    plt.tight_layout()
    if save is not None:
        if save is True and type(fname)==str:
            plt.savefig(fname+".png")
            plt.savefig(fname+".pdf")
        elif type(save)==str:
            plt.savefig(save+".png")
            plt.savefig(save+".pdf")
        else:
            print "Figure not saved, specify a name of the figure file"



