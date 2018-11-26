
# coding: utf-8

# In[1]:

text_size = 8
tick_size = 8
import numpy as np
import deepdish.io as dd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tlmModel import tlmModel
Es_normal={'rhs':"Ms",
        'n':(128,),
        'l':(64.0,),
        'bc':"periodic",
        'it':"pseudo_spectral",
        'dt':0.1,
        'analyze':False,
        'verbose':False,
        'setPDE':False}
m = tlmModel(Es=Es_normal,Ps='auto/tlm_set8.hdf5')
pts_per_inch=72.27       # this is a latex constant, don't change it.
text_width_in_pts=252.0  # write "\the\textwidth" (or "\showthe\columnwidth" for a 2 collumn text)
                         # inside a figure environment in latex, the result will be on the dvi/pdf next to the figure. See url above.
text_width_in_inches=text_width_in_pts/pts_per_inch
golden_ratio=0.618       # make rectangles with a nice proportion
#        squashed_ratio=0.309       # make rectangles with a nice proportion
inverse_latex_scale=2    # figure.png or figure.eps will be intentionally larger, because it is prettier
                         # when compiling latex code, use \includegraphics[scale=(1/inverse_latex_scale)]{figure}
fig_proportion = (2.0/3.0) # we want the figure to occupy 2/3 (for example) of the text width
csize=inverse_latex_scale*fig_proportion*text_width_in_inches
fig_size=(1.0*csize,golden_ratio*csize)  # always 1.0 on the first argument
text_size=inverse_latex_scale*text_size  # find out the fontsize of your latex text, and put it here
tick_size=inverse_latex_scale*tick_size
# learn how to configure: http://matplotlib.sourceforge.net/users/customizing.html
params = {#'backend': 'ps',
          'axes.labelsize': text_size,
          #'axes.linewidth' : 0,
          'font.size': text_size,
          'legend.fontsize': tick_size,
          'legend.handlelength': 2.5,
          'legend.borderaxespad': 0,
          'xtick.labelsize': tick_size,
          'ytick.labelsize': tick_size,
          'font.family':'serif',
          'font.size': text_size,
          #'font.serif':['Times'], # Times, Palatino, New Century Schoolbook, Bookman, Computer Modern Roman
          #'ps.usedistiller': 'xpdf',
          #'text.usetex': True,
          'figure.figsize': fig_size,
          #'text.latex.preamble' : [ r'\usepackage{amsmath}', # include here any neede package for latex
          #                        ] ,
          }
#        gs2 = gridspec.GridSpec(1, 1)
#        gs2.update(left=0.75, right=0.95, hspace=0.2)
#        ax3 = plt.subplot(gs2[0,0],adjustable='box-forced')
plt.rcParams.update(params)
figsize = (1.618*6,6)
get_ipython().magic(u'matplotlib inline')


# In[2]:

plt.figure(figsize=figsize)
print 100/m.p['conv_P']
t,chi02=m.ode_integrate([0.001,0.1,0.1],p=600/m.p['conv_P'],chi=0.2,beta=1.0,iota=0.0,finish=1*m.p['conv_T_to_t'])
t,chi05=m.ode_integrate([0.001,0.1,0.1],p=600/m.p['conv_P'],chi=0.5,beta=1.0,iota=0.0,finish=1*m.p['conv_T_to_t'])
t,chi07=m.ode_integrate([0.001,0.1,0.1],p=600/m.p['conv_P'],chi=0.7,beta=1.0,iota=0.0,finish=1*m.p['conv_T_to_t'])
plt.plot(t/m.p['conv_T_to_t'],chi02[0],'r',label=r'$\chi=0.2$')
plt.plot(t/m.p['conv_T_to_t'],chi05[0],'g',label=r'$\chi=0.5$')
plt.plot(t/m.p['conv_T_to_t'],chi07[0],'b',label=r'$\chi=0.7$')
plt.xlim([t[0]/m.p['conv_T_to_t'],t[-1]/m.p['conv_T_to_t']])
plt.xlabel(r'$T\,[yr]$')
plt.ylabel(r'$b$')
plt.legend(loc='best')
plt.tight_layout()
#plt.savefig('../../Dropbox/code/tlm/tlm_time_integration_chi_02_05_07.png')
#plt.savefig('../../Dropbox/code/tlm/tlm_time_integration_chi_02_05_07.pdf')

