#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
#  plot_final_field.py
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
#  
"""  
__version__= 1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

def read_col(filename, numcol, ndim):
    return np.loadtxt(filename, usecols = (numcol,)).reshape((ndim,ndim))

def plot_wood_herb_fields(filename):
    #ndim = 64
    data = np.loadtxt(filename, usecols = (0,))
    ndim = len(data)
    #print ndim
    ndim = np.sqrt(ndim)
    wood = read_col(filename, 0 , ndim)
    herb = read_col(filename, 1 , ndim)
    upper = read_col(filename, 2 , ndim)
    lower = read_col(filename, 3 , ndim)
    #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    ax1 = plt.subplot(2, 2, 1)
    cb1 = ax1.imshow(wood,cmap=plt.cm.YlGn)#, vmin=0.0, vmax=1.0)
    plt.colorbar(cb1)
    ax1.set_adjustable('box-forced')
    ax1.autoscale(False)
    ax1.set_title('Woody specie')
    ax2 = plt.subplot(2, 2, 2)
    cb2 = ax2.imshow(herb,cmap=plt.cm.YlGn)#, vmin=0.0, vmax=0.1)
    plt.colorbar(cb2)
    ax2.set_adjustable('box-forced')
    ax2.autoscale(False)
    ax2.set_title('Herb specie')
    ax3 = plt.subplot(2, 2, 4)
    cb3 = ax3.imshow(upper,cmap=plt.cm.Blues)#, vmin=6.2, vmax=6.3)
    plt.colorbar(cb3)
    ax3.set_adjustable('box-forced')
    ax3.autoscale(False)
    ax3.set_title('Upper layer')
    ax4 = plt.subplot(2, 2, 3)
    cb4 = ax4.imshow(lower,cmap=plt.cm.Blues)#, vmin=13.8, vmax=13.9)
    plt.colorbar(cb4)
    ax4.set_adjustable('box-forced')
    ax4.autoscale(False)
    ax4.set_title('Lower layer')
    
    #plt.savefig(filename+".png", fmt="png")
    plt.show()


def plot_one_layer_fields(filename):
    data = np.loadtxt(filename, usecols = (0,))
    ndim = len(data)
    #print ndim
    ndim = np.sqrt(ndim)
    wood = read_col(filename, 0 , ndim)
    upper = read_col(filename, 1 , ndim)
    lower = read_col(filename, 3 , ndim)
    ax1 = plt.subplot(2,2,1)
    #f, (ax1, ax2,ax3) = plt.subplots(3, sharex='col', sharey='row')
    cb1 = ax1.imshow(wood,cmap=plt.cm.YlGn)#, vmin=0, vmax=1)
    plt.colorbar(cb1)
    ax1.set_adjustable('box-forced')
    ax1.autoscale(False)
    ax1.set_title('Woody specie')
    ax2 = plt.subplot(2,2,2)
    cb2 = ax2.imshow(upper,cmap=plt.cm.YlGn)#, vmin=6.2, vmax=6.3)
    plt.colorbar(cb2)
    ax2.set_adjustable('box-forced')
    ax2.autoscale(False)
    ax2.set_title('Herb specie')
    ax3 = plt.subplot(2,2,4)
    cb3 = ax3.imshow(lower,cmap=plt.cm.Blues)#, vmin=13.8, vmax=13.9)
    plt.colorbar(cb3)
    ax3.set_adjustable('box-forced')
    ax3.autoscale(False)
    ax3.set_title('One layer')
    plt.tight_layout()
    plt.show()

def plot_wood_fields(filename):
    data = np.loadtxt(filename, usecols = (0,))
    ndim = len(data)
    #print ndim
    ndim = np.sqrt(ndim)
    wood = read_col(filename, 0 , ndim)
    upper = read_col(filename, 2 , ndim)
    lower = read_col(filename, 3 , ndim)
    ax1 = plt.subplot(2,2,1)
    #f, (ax1, ax2,ax3) = plt.subplots(3, sharex='col', sharey='row')
    cb1 = ax1.imshow(wood,cmap=plt.cm.YlGn)#, vmin=0, vmax=1)
    plt.colorbar(cb1)
    ax1.set_adjustable('box-forced')
    ax1.autoscale(False)
    ax1.set_title('Woody specie')
    ax2 = plt.subplot(2,2,3)
    cb2 = ax2.imshow(upper,cmap=plt.cm.Blues)#, vmin=0, vmax=10)
    plt.colorbar(cb2)
    ax2.set_adjustable('box-forced')
    ax2.autoscale(False)
    ax2.set_title('Upper layer')
    ax3 = plt.subplot(2,2,4)
    cb3 = ax3.imshow(lower,cmap=plt.cm.Blues)#, vmin=0, vmax=30)
    plt.colorbar(cb3)
    ax3.set_adjustable('box-forced')
    ax3.autoscale(False)
    ax3.set_title('Lower layer')
    plt.tight_layout()
    plt.show()

def plot_herb_fields(filename):
    data = np.loadtxt(filename, usecols = (0,))
    ndim = len(data)
    #print ndim
    ndim = np.sqrt(ndim)
    wood = read_col(filename, 1 , ndim)
    upper = read_col(filename, 2 , ndim)
    lower = read_col(filename, 3 , ndim)
    ax1 = plt.subplot(2,2,2)
    #f, (ax1, ax2,ax3) = plt.subplots(3, sharex='col', sharey='row')
    cb1 = ax1.imshow(wood,cmap=plt.cm.YlGn)#, vmin=0, vmax=1)
    plt.colorbar(cb1)
    ax1.set_adjustable('box-forced')
    ax1.autoscale(False)
    ax1.set_title('Herb specie')
    ax2 = plt.subplot(2,2,3)
    cb2 = ax2.imshow(upper,cmap=plt.cm.Blues)#, vmin=0, vmax=10)
    plt.colorbar(cb2)
    ax2.set_adjustable('box-forced')
    ax2.autoscale(False)
    ax2.set_title('Upper layer')
    ax3 = plt.subplot(2,2,4)
    cb3 = ax3.imshow(lower,cmap=plt.cm.Blues)#, vmin=0, vmax=30)
    plt.colorbar(cb3)
    ax3.set_adjustable('box-forced')
    ax3.autoscale(False)
    ax3.set_title('Lower layer')
    plt.tight_layout()
    plt.show()

        
def plot_vegetation(filename):
    data = np.loadtxt(filename, usecols = (0,))
    ndim = len(data)
    ndim = np.sqrt(ndim)
    wood = read_col(filename, 0 , ndim)
    herb = read_col(filename, 1 , ndim)    
    fig, (ax,ax2) = plt.subplots(1,2,sharey=True)
    ax.tick_params(axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    labelbottom='off')
    ax2.tick_params(axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    labelbottom='off')
    ax.tick_params(axis='y',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    left='on',      # ticks along the bottom edge are off
                    right='off',         # ticks along the top edge are off
                    labelbottom='off')
    ax2.tick_params(axis='y',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    left='off',      # ticks along the bottom edge are off
                    right='off',         # ticks along the top edge are off
                    labelbottom='off')
    #ax = fig.add_subplot(1, 2, 1)
    ax.imshow(wood,cmap=plt.cm.YlGn, vmin=0.0, vmax=0.5)
    ax.set_adjustable('box-forced')
    #ax.autoscale(False)
    #fig.colorbar(ax, orientation='horizontal')
    #cb1 = plt.colorbar(ax,shrink=0.5)
    #plt.title('Woody specie')
    #ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_adjustable('box-forced')
    ax2.imshow(herb,cmap=plt.cm.YlGn, vmin=0.0, vmax=0.05)
    #ax2.autoscale(False)
    #fig.colorbar(ax2, orientation='horizontal')
    #cb2 = plt.colorbar(ax2,shrink=0.5)
    #plt.title('Herb specie')
    plt.savefig(filename+".png",bbox_inches='tight')
    #plt.show()
    
def plot_water(upper,lower):
    fig = plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    ax.imshow(upper,cmap=plt.cm.Blues)#, vmin=0.0, vmax=0.5)
    ax.set_adjustable('box-forced')
    ax.autoscale(False)
    #cb1 = plt.colorbar(ax,shrink=0.5)
    plt.title('Upper layer')
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.set_adjustable('box-forced')
    ax2.imshow(lower,cmap=plt.cm.Blues)#, vmin=0.0, vmax=1.0)
    ax2.autoscale(False)
    #cb2 = plt.colorbar(ax2,shrink=0.5)
    plt.title('Lower layer')
    plt.savefig("waters.pdf", format='pdf')
    #plt.show()
    

def main(args):
    if args.b2wh:
        plot_one_layer_fields(args.filename)
    elif args.btw2h:
        plot_wood_fields(args.filename)
    elif args.bhw2h:
        plot_herb_fields(args.filename)
    elif args.one_layer:
        plot_one_layer_fields(args.filename)
    elif args.veg:
        plot_vegetation(args.filename)
    else:
        plot_wood_herb_fields(args.filename)
        
    


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

