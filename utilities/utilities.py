# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:04:05 2016

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""
import numpy as np
import scipy.spatial as spat

def ReadSimulationResult(fname,file_format="npy"):
    if file_format == "npy":
        data=np.load(fname+".npy")
    elif file_format == "txt":
        data=np.loadtxt(fname+".dat",skiprows=1 )
    else:
        print "Uknown format"
    return data
    
def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx
    
def sort_dots(x,y, start=0,metrics="mahalanobis"):
    data = np.array([x,y])
    dist_m = spat.distance.squareform(spat.distance.pdist(data.T, metrics))

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