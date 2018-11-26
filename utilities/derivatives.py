# -*- coding: utf-8 -*-
"""
# Copyright (C) 2016 by
#    Omer Tzuk    <omertz@post.bgu.ac.il>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY <copyright holder> ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import numpy as np
from numba import guvectorize, float64
__version__=1.0
__author__ = """Omer Tzuk (omertz@post.bgu.ac.il)"""
from numba.pycc import CC
#
cc = CC('derivatives_comp')
## Uncomment the following line to print out the compilation steps
cc.verbose = True

""" Laplacian functions """
""" 1D                  """
@cc.export('neumann_laplacian_1d', 'f8[:](f8[:], f8)')
def neumann_laplacian_1d(u,dx2):
    """Return finite difference Laplacian approximation of 1d array.
    Uses Neumann boundary conditions and a 2nd order approximation.
    """
    nx, = u.shape
    aa = np.zeros(nx+2)
    aa[0]=u[1]
    aa[-1]=u[-2]
    aa[1:nx+1] = u
    laplacian=(aa[0:nx]+aa[2:nx+2]-(2.0)*aa[1:nx+1])/dx2
    return laplacian



@guvectorize([(float64[:], float64, float64[:])], '(n),()->(n)')
def neumann_laplacian_1d_jit(u,dx2,laplacian):
    """Return finite difference Laplacian approximation of 1d array.
    Uses Neumann boundary conditions and a 2nd order approximation.
    """
    nx, = u.shape
    aa = np.zeros(nx+2)
    aa[0]=u[1]
    aa[-1]=u[-2]
    aa[1:nx+1] = u
    laplacian=(aa[0:nx]+aa[2:nx+2]-(2.0)*aa[1:nx+1])/dx2
#    return laplacian

#@autojit(target="cpu")
def periodic_laplacian_1d(u,dx2):
    """Return finite difference Laplacian approximation of 1d array.
    Uses periodic boundary conditions and a 2nd order approximation.
    """
    nx, = u.shape
    aa = np.zeros(nx+2)
    aa[0]=u[-1]
    aa[-1]=u[0]
    aa[1:nx+1] = u
    laplacian=(aa[0:nx]+aa[2:nx+2]-(2.0)*aa[1:nx+1])/dx2
    return laplacian

""" 2D                  """

@guvectorize([(float64[:,:], float64, float64[:,:])], '(n,n),()->(n,n)')
def neumann_laplacian_2d_jit(u,dx2,laplacian):
    """Return finite difference Laplacian approximation of 2d array.
    Uses Neumann boundary conditions and a 2nd order approximation.
    """
    nx, ny = u.shape
    aa = np.zeros((nx+2,ny+2))
    aa[1:nx+1,0]=u[:,1]
    aa[1:nx+1,-1]=u[:,-2]
    aa[0,1:ny+1]=u[1,:]
    aa[-1,1:ny+1]=u[-2,:]
    aa[1:nx+1,1:ny+1] = u
    laplacian=(aa[0:nx,1:ny+1]+aa[2:nx+2,1:ny+1]+aa[1:nx+1,0:ny]+aa[1:nx+1,2:ny+2]-(4.0)*aa[1:nx+1,1:ny+1])/dx2
#    return laplacian
#@cc.export('neumann_laplacian_2d', 'f8[:,:](f8[:,:], f8)')
def neumann_laplacian_2d(u,dx2):
    """Return finite difference Laplacian approximation of 2d array.
    Uses Neumann boundary conditions and a 2nd order approximation.
    """
    nx, ny = u.shape
    aa = np.zeros((nx+2,ny+2))
    aa[1:nx+1,0]=u[:,1]
    aa[1:nx+1,-1]=u[:,-2]
    aa[0,1:ny+1]=u[1,:]
    aa[-1,1:ny+1]=u[-2,:]
    aa[1:nx+1,1:ny+1] = u
    laplacian=(aa[0:nx,1:ny+1]+aa[2:nx+2,1:ny+1]+aa[1:nx+1,0:ny]+aa[1:nx+1,2:ny+2]-(4.0)*aa[1:nx+1,1:ny+1])/dx2
    return laplacian

#@autojit(target="cpu")
def periodic_laplacian_2d(u,dx2):
    """Return finite difference Laplacian approximation of 2d array.
    Uses periodic boundary conditions and a 2nd order approximation.
    """
    nx, ny = u.shape
    aa = np.zeros((nx+2,ny+2))
    aa[1:nx+1,0]=u[:,-1]
    aa[1:nx+1,-1]=u[:,0]
    aa[0,1:ny+1]=u[-1,:]
    aa[-1,1:ny+1]=u[0,:]
    aa[1:nx+1,1:ny+1] = u
    laplacian=(aa[0:nx,1:ny+1]+aa[2:nx+2,1:ny+1]+aa[1:nx+1,0:ny]+aa[1:nx+1,2:ny+2]-(4.0)*aa[1:nx+1,1:ny+1])/dx2
    return laplacian

""" Gradient functions """

#@autojit(target="cpu")
def neumann_gradient_1d(array1D, dx=1.0):
	""" (1D array, dx^2) -> d(1D array)/dx
	Implementing the 2nd order 1D partial_x derivative with neumann condition
	"""
	gradient = np.zeros_like(array1D)
	gradient[1:-1] = ((1.0/2.0)*array1D[2:]-(1.0/2.0)*array1D[:-2])/ dx
	gradient[0] = 0.0
	gradient[-1]= 0.0
	return gradient

#@autojit(target="cpu")
def periodic_gradient_1d(u,dx=1):
    """Return finite difference Gradient approximation of 2d array.

    Uses periodic boundary conditions and a 2nd order approximation.

    """
    gradient = (  np.roll(u,-1,axis=0)
                 - np.roll(u, 1,axis=0) )/(2*dx)

    return gradient

#@autojit(target="cpu")
def neumann_gradient_2d(u,dx=1):
    """Return finite difference Gradient approximation of 2d array.

    Uses Neumann boundary conditions and a 2nd order approximation.

    """
    n=len(u[:,0])
    gradient = np.zeros([n,n+2])
    gradient[:,1:-1] = u
    gradient[:,0] = u[:,1]
    gradient[:,-1] = u[:,-2]
    gradient[:,1:-1] = (  gradient[:,:-2]
                          - gradient[:,2:] )/(2*dx)

    return gradient[:,1:-1]

#@autojit(target="cpu")
def periodic_gradient_2d(u,dx=1):
    """Return finite difference Gradient approximation of 2d array.

    Uses Neumann boundary conditions and a 2nd order approximation.

    """
    n=len(u[:,0])
    gradient = np.zeros([n,n+2])
    gradient[:,1:-1] = u
    gradient[:,0] = u[:,1]
    gradient[:,-1] = u[:,-2]
    gradient[:,1:-1] = (  gradient[:,:-2]
                          - gradient[:,2:] )/(2*dx)

    return gradient[:,1:-1]



if __name__ == "__main__":
#    pass
    cc.compile()