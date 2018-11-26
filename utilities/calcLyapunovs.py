# 
# See the exercise 
# http://www.physics.cornell.edu/~myers/teaching/ComputationalMethods/ComputerExercises/ChaosLyapunov.html
#
# <nbformat>2</nbformat>
# <codecell>
"""Chaos and Lyapunov exercise"""

# <codecell>
from IterateLogistic import *
import scipy.optimize

# <markdowncell>
# Calculate Lyapunov exponent for mu=0.9
#
# traj = TrajectoryDifference(f,0.1,0.1000000000001, 100, 0.9)
# p = FitLyapunovExponent(traj)
# p = p[0] # Weird error message removal
# PlotFit(traj, p)

# <codecell>
def TrajectoryDifference(g, x1, x2, N, mu):
    """
    Calculates the difference between two trajectories that start
    at two points x1 and x2, presumably close to one another.
    Returns list of length N+1,
    [x2-x1, g(x2)-g(x1), g(g(x2))-g(g(x1)), ...]

    Used for illustrating sensitive dependence on initial conditions for
    chaotic regions, and for calculating Lyapunov exponents.

    You can usefully pick x2-x1 comparable to machine precision (seventeen
    digits). Don't choose x1 or x2 to equal zero or one in the logistic 
    map: zero is an unstable fixed point.
    """
    pass

# <codecell>
def PlotTrajectoryDifference(g, x1, x2, N, mu):
    """
    Calls TrajectoryDifference to find the difference, then plots the 
    absolute value of the difference using 
        pylab.semilogy(scipy.fabs(dx)).
    (Given just one array, pylab assumes the other axis is just the
    index into the array.)

    Notice that the differences stop growing when they become of order 
    one (naturally). Don't use such long trajectories to calculate the 
    Lyapunov exponents: it will bias the results downward.
    """
    pass

# <codecell>
def LyapunovFitFunc(p, traj_diff):
    """
    Given a trajectory difference traj_diff and a tuple
        p = (lyapExponent, lyapLogPrefactor)
    returns the residuals 
        log(|y_n|) - log( exp(lyapLogPrefactor + lyapExponent * n) )
     == log(|y_n|) - (lyapLogPrefactor + lyapExponent * n) (avoids overflow)

    The residual is the difference between the data and the fit:
    if the residuals were zero, the trajectory difference would be 
    perfectly described as a growing exponential.
    That is, the growth of the difference between two trajectories is 
    expected to be of the form
        |x_n - y_n| \sim lyapPrefactor * exp(lyapExponent n)

    We take the log of the difference so that the least--squares fit
    will emphasize the initial points and final points roughly equally.
    We use the log of the lyapPrefactor because it must be positive:
    it's a standard trick in nonlinear fitting to change variables like
    this to enforce ranges in parameters.

    Used by FitLyapunovExponent to generate a least-squares fit for 
    the Lyapunov exponent and prefactor.
    """
    pass

# <codecell>
def FitLyapunovExponent(traj_diff, p0=(1., -13.)):
    """
    Given a trajectory difference and an initial guess 
    p0=(lyapExponent, lyapPrefactor) for the Lyapunov exponent and 
    prefactor, uses scipy.optimize.leastsq to do a best fit for 
    these two constants, and returns them.

    The return value will likely include a warning message, even though
    the fit seems fine and the warning meaningless. You'll likely need 
    to delete it before using PlotFit.
    """
    pass

# <codecell>
def PlotFit(traj_diff, p):
    """
    Given a trajectory difference and p=(lyapExponent, lyapPrefactor),
    plot |traj_diff| and the fit on a semilog y axis.
    """
    pass


# <markdowncell>
# Copyright (C) Cornell University
# All rights reserved.
# Apache License, Version 2.0


