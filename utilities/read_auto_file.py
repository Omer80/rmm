# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:29:58 2016

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""
import numpy as np
import re
from astropy.table import Table
import matplotlib.pyplot as plt
import scipy.spatial as ss
import itertools
#import matplotlib.cm as cm
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
#table = sort_bifurcation(read_b("./auto/b.b2s1_data"), "b2s1_ode")
def readBf_b(fname):
    lines = (line for line in open(fname, 'rb')  if re.match(r'^\s*[1-9]', line) )
    data = np.genfromtxt(lines,unpack=True)
    return data

def sortBf(data, model):
    if model=="b2s1_ode":
        col_names = ['branch','stability','p','L2','b1','b2','s']
#        stability = data[1]>=0
        data_cols = [data[0],data[1],data[4],data[5],data[6],data[7],data[8]]
        t = Table(data_cols, names=col_names
                  ,meta={'name': model+' ODE AUTO bifurcation diagram'}
                  ,dtype=('i4', 'i4','f8','f8','f8','f8','f8'))
        t.add_index('branch')
    elif model=="b2s2":
        col_names = ['branch','stability','p','L2','b1','b2','s1','s2']
#        stability = data[1]>=0
        data_cols = [data[0],data[1],data[4],data[5],data[6],data[7],data[8],data[9]]
        t = Table(data_cols, names=col_names
                  ,meta={'name': model+' ODE AUTO bifurcation diagram'}
                  ,dtype=('i4', 'i4','f8','f8','f8','f8','f8','f8'))
        t.add_index('branch')
    return t

def save_h5(table,fname,pather,overwriter=True):
    table.write(fname+".hdf5", path=pather, append=True,overwrite=overwriter)

def load_h5(fname, pather):
    table = Table.read(fname+".hdf5", path=pather)
    return table

def plotBf(table, control_parameter, states, n_branch):
    f, ax = plt.subplots(len(states), sharex=True)
    ax[0].set_title(r'Bifurcation diagram')
    table.add_index("branch")
    for i,state in enumerate(states):
        colors = itertools.cycle(["r", "b", "g", "m", "c"])
        ax[i].set_ylabel(str(state))
        for branch in n_branch:
            cl = next(colors)
            x = table.loc[branch][control_parameter]
            y = table.loc[branch][state]
            sorted_index = sort_dots(x,y)
            x,y = radial_sort_line(x[sorted_index],y[sorted_index])
    #        x,y = radial_sort_line(x,y)
            stable   = np.greater(table.loc[branch]['stability'],0)
            unstable = np.less(table.loc[branch]['stability'],0)
    #        plt.plot(x[sorted_index],y[sorted_index],color=cl)
            ax[i].plot(np.ma.masked_where(unstable, x),np.ma.masked_where(unstable, y), color=cl,linestyle="-")
        plt.plot(np.ma.masked_where(stable, x),np.ma.masked_where(stable, y), color=cl,linestyle="--")
    ax[-1].set_xlabel(r'p')
    plt.show()
#    return fig

def distance(pt_1, pt_2):
    pt_1 = np.array((pt_1[0], pt_1[1]))
    pt_2 = np.array((pt_2[0], pt_2[1]))
    return np.linalg.norm(pt_1-pt_2)

def radial_sort_line(x, y):
    """Sort unordered verts of an unclosed line by angle from their center."""
    # Radial sort
    x0, y0 = x.mean(), y.mean()
    angle = np.arctan2(y - y0, x - x0)

    idx = angle.argsort()
    x, y = x[idx], y[idx]

    # Split at opening in line
    dx = np.diff(np.append(x, x[-1]))
    dy = np.diff(np.append(y, y[-1]))
    max_gap = np.abs(np.hypot(dx, dy)).argmax() + 1

    x = np.append(x[max_gap:], x[:max_gap])
    y = np.append(y[max_gap:], y[:max_gap])
    return x, y

def sort_dots(x,y, start=0,metrics="mahalanobis"):
    data = np.array([x,y])
    dist_m = ss.distance.squareform(ss.distance.pdist(data.T, metrics))

    total_points = data.shape[1]
    points_index = set(range(total_points))
    sorted_index = []
    target    = start

    points_index.discard(target)
    while len(points_index)>0:
        candidate = list(points_index)
        nneigbour = candidate[dist_m[target, candidate].argmin()]
        points_index.discard(nneigbour)
        points_index.discard(target)
        #print points_index, target, nneigbour
        sorted_index.append(target)
        target    = nneigbour
    sorted_index.append(target)
    return sorted_index

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)