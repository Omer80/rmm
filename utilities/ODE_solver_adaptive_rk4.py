#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 19:03:38 2017

@author: ohm
"""

#    Adaptive step size ODE solver with RK5 (embedded RK4)
#
#    Formulation and notation taken from Numerical Recipes in Fortran 
#    77, Second Edition (1992) by William H. Press, Brian P. Flannery, 
#    Saul A. Teukolsky, and William T. Vetterling.
#    ISBN-10: 052143064X 
#
#    Alexander Miles - Last edit: 7/30/12

from time import time # Import time-keeping library
tstart = time()       # Define starting time
print "[--.--] Starting. Importing libraries."
from pylab import *   # Import libraries for plotting results
close('all')          # Close previously opened plots


### Cash-Karp Parameters - From literature

print "[%4.3f] Defining parameters and functions." % (time()-tstart)

a2,   a3,  a4,  a5,  a6      =        1/5.,    3/10.,       3/5.,            1.,        7/8.
b21, b31, b32, b41, b42, b43 =        1/5.,    3/40.,      9/40.,         3/10.,      -9/10., 6/5.
b51, b52, b53, b54           =     -11/54.,     5/2.,    -70/27.,        35/27.
b61, b62, b63, b64, b65      = 1631/55296., 175/512., 575/13824., 44275/110592.,   253/4096. 
c1,   c2,  c3,  c4,  c5, c6  =     37/378.,       0.,   250/621.,      125/594.,          0.,  512/1771.
c1star, c2star, c3star, c4star, c5star, c6star = 2825/27648., 0.,  18575/48384.,13525/55296., 277/14336., 1/4.

def control(deriv, tol, y0, t0, h, tmax, v=False):
    '''
    This funciton takes in a python function that returns the derivative,
    a tolerance for error between RK5 and RK4, initial conditions on y (dependant) 
    and t (independant) as well as an initial step size.
    
    Keyword arguments:
    v - Verbose
    '''
    tstart = time()
    if v==True: print "[%4.3f] Solving with initial condition (%0.2f, %0.2f), step size of %0.2f from t=0...%0.2f" % ((time()-tstart), y0, t0, h, tmax)

    y = [y0]      # Set up the initial conditions on y
    t = [t0]      # and t while creating the output lists
   
    t_curr, y_curr, count, ncount = t0, y0, 0, 0 # Setup counters and trackers

    while t_curr < tmax:
        t_next, y_next, h, h1 = stepper(deriv, t_curr, y_curr, h, tol)
        if h1 < 0.9*h: 
            if v==True: print "[%4.3f] Reduced step size from %0.4e to %0.4e at t = %0.2f" % ((time()-tstart),h, h1, t_curr)
            h = h1
        elif h1 > 1.1*h:
            if v==True: print "[%4.3f] Increased step size from %0.4e to %0.4e at t = %0.2f" % ((time()-tstart),h, h1, t_curr)
            h = h1
        else:
            y.append(y_next)
            t.append(t_next)
            y_curr, t_curr = y_next, t_next
            ncount += 1
        count += 1
    if v==True: print "[%4.3f] Done. %i iterations, %i points" % ( (time()-tstart), count, ncount )
    return y, t

def stepper(deriv, t, y, h, tol):
   '''
   This function is called by the control function to take
   a single step forward. The inputs are the derivative function,
   the previous time and function value, the step size in time (h),
   and the tolerance for error between 5th order Runge-Kutta and 4th
   order Runge-Kutta.
   '''

   k1 = h*deriv(t,y)
   k2 = h*deriv(t+a2*h,y+b21*k1)
   k3 = h*deriv(t+a3*h,y+b31*k1+b32*k2)
   k4 = h*deriv(t+a4*h,y+b41*k1+b42*k2+b43*k3)
   k5 = h*deriv(t+a5*h,y+b51*k1+b52*k2+b53*k3+b54*k4)
   k6 = h*deriv(t+a6*h,y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5)
   y_n_plus_1      = y +     c1*k1 +     c2*k2 +     c3*k3 +     c4*k4 +     c5*k5 +     c6*k6
   y_n_plus_1_star = y + c1star*k1 + c2star*k2 + c3star*k3 + c4star*k4 + c5star*k5 + c6star*k6
   DELTA           = y_n_plus_1 - y_n_plus_1_star
   try:
       h1 = h*abs(tol/DELTA)**0.2    # Finds step size required to meet given tolerance
   except ZeroDivisionError:
       h1 = h                        # When you are very close to ideal step, DELTA can be zero
   return t+h, y_n_plus_1, h, h1


#### From here on down is a test case. Comment out for actual use.

def dydt(t,y):
    gamma = 2.0
    return -gamma*y

y0   = 1.0           # Initial conditions
t0   = 0.0           #
h    = 0.1           # Initial step size
tmax = 2.0           # End of time interval
npts = int(tmax/h)   # Number of points
tol  = 1E-12         # The desired error between 4th and 5th order
                     # note that this IS NOT the error between numeric
                     # solution and the actual solution

a, b = control(dydt, tol, y0, t0, h, tmax, v=True)
plot(b,a)
show()