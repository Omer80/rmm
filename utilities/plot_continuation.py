# -*- coding: utf-8 -*-
"""
#  b2s2_fortran_continuation.py
#  
#  Copyright 2015 Omer Tzuk <omertz@post.bgu.ac.il>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
"""
__version__= 1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
import argparse
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def readhdf(hdf_fname):
    """ Get the hdf5 filename, and returns the continuation data"""
    f = h5py.File("%s.hdf5"%hdf_fname,'r')
    cont = f["continuation"]
    prec = cont["prec"][:]
    cont_b1 = cont["b1"]
    cont_b2 = cont["b2"]
    cont_s1 = cont["s1"]
    cont_s2 = cont["s2"]
    integral_b1 = cont_b1["integral"][:]
    l2norm_b1 = cont_b1["l2norm"][:]
    max_b1 = cont_b1["max"][:]
    mean_b1 = cont_b1["mean"][:]
    min_b1 = cont_b1["min"][:]
    b1 = {"integral":integral_b1, "l2norm":l2norm_b1, "max":max_b1, "mean":mean_b1, "min":min_b1}
    integral_b2 = cont_b2["integral"][:]
    l2norm_b2 = cont_b2["l2norm"][:]
    max_b2 = cont_b2["max"][:]
    mean_b2 = cont_b2["mean"][:]
    min_b2 = cont_b2["min"][:]
    b2 = {"integral":integral_b2, "l2norm":l2norm_b2, "max":max_b2, "mean":mean_b2, "min":min_b2}
    integral_s1 = cont_s1["integral"][:]
    l2norm_s1 = cont_s1["l2norm"][:]
    max_s1 = cont_s1["max"][:]
    mean_s1 = cont_s1["mean"][:]
    min_s1 = cont_s1["min"][:]
    s1 = {"integral":integral_s1, "l2norm":l2norm_s1, "max":max_s1, "mean":mean_s1, "min":min_s1}
    integral_s2 = cont_s2["integral"][:]
    l2norm_s2 = cont_s2["l2norm"][:]
    max_s2 = cont_s2["max"][:]
    mean_s2 = cont_s2["mean"][:]
    min_s2 = cont_s2["min"][:]
    s2 = {"integral":integral_s2, "l2norm":l2norm_s2, "max":max_s2, "mean":mean_s2, "min":min_s2}
    data = {"prec":prec, "b1":b1, "b2":b2, "s1":s1, "s2":s2}
    return data

def plot_continuation(data):
    """receive """
    f, axarr = plt.subplots(4, sharex=True)
    P = data['prec']
    axarr[0].plot(P,data['b1']['mean'], 'g.',P,data['b1']['min'], 'g--',P,data['b1']['max'], 'g-')
    axarr[0].set_title('Biomass and ground water saturation')
    axarr[0].set_ylabel(r"$b_1$")
    axarr[1].plot(P,data['b2']['mean'],'g.',P,data['b2']['min'],'g--',P,data['b2']['max'],'g-')
    axarr[1].set_ylabel(r"$b_2$")
    axarr[2].plot(P,data['s1']['mean'], 'b.',P,data['s1']['min'], 'b--',P,data['s1']['max'], 'b-')
    axarr[2].set_ylabel(r"$s_1$")
    axarr[3].plot(P,data['s2']['mean'],'b.',P,data['s2']['min'],'b--',P,data['s2']['max'],'b-')
    axarr[3].set_ylabel(r"$s_2$")
    axarr[3].set_xlabel(r"$prec$")
    plt.show()

def main(args):
    print "Everything's sweet!"
    data = readhdf("continuation")
    plot_continuation(data)
    return 0

def add_parser_arguments(parser):
    "Add command-line arguments."""
    parser.add_argument('filename', type=str, nargs='?', default=None,
                        help='output file')
    parser.add_argument("--b2wh",
                        action="store_true",
                        dest="b2wh",
                        default=False,
                        help="Plot one layer fields")
    parser.add_argument("--btw2h",
                        action="store_true",
                        dest="btw2h",
                        default=False,
                        help="Plot wood and water fields")
    parser.add_argument("--veg",
                        action="store_true",
                        dest="veg",
                        default=False,
                        help="Plot vegetation fields")
    parser.add_argument("--bhw2h",
                        action="store_true",
                        dest="bhw2h",
                        default=False,
                        help="Plot wood and water fields")
    parser.add_argument("--one_layer",
                        action="store_true",
                        dest="one_layer",
                        default=False,
                        help="Plot one layer model fields")
    return parser
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    args = parser.parse_args()
    main(args)