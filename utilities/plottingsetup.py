#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 18:26:14 2017

@author: ohm
"""

import numpy as np
import matplotlib as mpl
#mpl.use('pgf')
#import itertools

def figsize(scale):
    #fig_width_pt = 345.0                          # Get this from LaTeX using \the\textwidth
    fig_width_pt = 379.4175                          # PloS Journal
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size
# learn how to configure: http://matplotlib.sourceforge.net/users/customizing.html
pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "font.size": 10,
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)
import matplotlib.pyplot as plt

def newfig(width):
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    return fig, ax


def plot_b_to_chi_per_p(p,iota,beta=1.0,initial_state=[0.9,0.2,0.2],
                  legend=None,save=None):
    from scanGrids import scanChigradient
    fig,ax=newfig(1.0)
    if type(p)!=list:
        p = [p]
    cl = plt.cm.winter(np.linspace(1, 0, len(p)+2))[1:-1]
    for i,p_i in enumerate(p):
        data = scanChigradient(p=p_i,iota=iota,beta=beta,initial_conditions=initial_state)
        plt.plot(data['chi'],data['mean'],color=cl[i], lw=2,label=r'$p={:3.2f}$'.format(p_i))
    plt.xlabel(r'$\chi$')
    plt.ylabel(r'$\left<b\right>$')
    #plt.ylim([-0.05,1.0])
    plt.xlim([0,1.0])
    #plt.axhline(0,color='k',linestyle='--',linewidth=1)
    if legend is not None:
        plt.legend(loc=legend,fontsize=13)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save+'.eps')
        #plt.savefig(save+'.png')


def plot_b_to_chi_per_a(p,a,beta=1.0,initial_state=[0.9,0.2,0.2],
                  legend=True,save=None):
    from scanGrids import scanChigradient
    fig,ax=newfig(1.0)
    if type(a)!=list:
        a = [a]
    cl = plt.cm.gist_heat(np.linspace(1, 0, len(a)+2))[1:-1]
    for i,a_i in enumerate(a):
        data = scanChigradient(p=p,iota=a_i,beta=beta,initial_conditions=initial_state)
        plt.plot(data['chi'],data['mean'],color=cl[i], lw=2,label=r'$a={:3.2f}$'.format(a_i))
    plt.xlabel(r'$\chi$')
    plt.ylabel(r'$\left<b\right>$')
    #plt.ylim([-0.05,1.0])
    plt.xlim([0,1.0])
    #plt.axhline(0,color='k',linestyle='--',linewidth=1)
    if legend is not None:
        plt.legend(loc=legend,fontsize=13)
    plt.tight_layout()
    if save is not None:
        plt.savefig(save+'.eps')
        #plt.savefig(save+'.png')