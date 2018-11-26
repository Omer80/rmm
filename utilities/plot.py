#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 14:32:27 2017

@author: ohm
"""
text_size = 12
tick_size = 12
import matplotlib.pyplot as plt
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