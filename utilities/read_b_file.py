# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:29:58 2016

@author: Omer Tzuk <omertz@post.bgu.ac.il>
"""

def read_b(fname):
    with open(fname) as f:
        for line in f:
            print line[0]