# -*- coding: utf-8 -*-
"""
#  rmModel.py
#
#  Rosenzweig - MacArthur Model
#
#  Copyright 2016 Omer Tzuk <omertz@post.bgu.ac.il>
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
import time
from sympy import symbols, Matrix,lambdify
from sympy.utilities.autowrap import ufuncify
import numpy as np
import scipy.linalg as linalg
from scipy.integrate import solve_ivp
from scipy.optimize import newton_krylov
from scipy.optimize.nonlin import NoConvergence
from scipy.optimize import root as root_ode
from scipy.fftpack import fftn, ifftn
from scipy.stats import norm,chi,pearson3
import scipy.sparse as sparse
from utilities import handle_netcdf_py3 as hn
import deepdish.io as dd

Es_normal={'rhs':"RM_Allee_effect",
        'n':(1024,),
        'l':(256.0,),
        'bc':"neumann",
        #'stochastic':("norm_around_p",0.03,100),
        'it':"scipy",
        'dt':0.001,
        'analyze':False,
        'verbose':True,
        'setPDE':False}

def main():
    global m,p
    m = rmModel(Es=Es_normal,Ps='auto/rm_set1.hdf5',Vs=None)
    return 0

class rmModel(object):
    def __init__(self,Ps,Es,Vs=None):
        if type(Ps)==str:
            self.Psfname=Ps
            self.p=dd.load(Ps)
        else:
            self.Psfname=None
            self.p = Ps
        self.setup=Es
        self.setup['nvar']=2
#        self.Vs=Vs
        self.verbose=Es['verbose']
        if self.verbose:
            start=time.time()
        self.set_equations()
        self.dt = Es['dt']
        self.time_elapsed = 0
        if self.setup['setPDE']:
            self.rhs=self.rhs_pde        
            self.p['nd']=len(Es['n'])
            if self.p['nd']==2:
                self.p['nx'],self.p['ny']=Es['n']
                self.p['lx'],self.p['ly']=Es['l']
                self.l=[self.p['lx'],self.p['ly']]
                self.n=[self.p['nx'],self.p['ny']]
                self.dg  = tuple([l/float(n) for l,n in zip(self.l,self.n)])
                self.dx  = self.dg[0]
            elif self.p['nd']==1:
                self.dg=[Es['l'][0]/float(Es['n'][0])]
                self.dx=self.dg[0]
            self.dx2 = self.dx**2
            self.dt=Es['dt']*self.dx2 / self.p['delta_s']
            self.X = np.linspace(0,Es['l'][0],Es['n'][0])
            from utilities.laplacian_sparse import create_laplacian #,create_gradient
            self.lapmat=create_laplacian(self.setup['n'],self.setup['l'], self.setup['bc'] , [1.0,self.p['delta_s'],self.p['delta_s']],verbose=self.verbose)
#            self.gradmat=create_gradient(self.setup['n'],self.setup['l'], self.setup['bc'] , [1.0,self.p['Dw'],self.p['Dh']])
            if self.verbose:
                print("Laplacian created")
        else:
            self.rhs=self.rhs_ode
        self.set_integrator()
        if Vs is not None:
            self.setup_initial_condition(Vs)
        if self.verbose:
            print("Time to setup: ",time.time()-start)
    """ Setting up model equations """
    def set_equations(self):
        x,y,t = symbols('x y t')
        self.Ps_symbols={}
        for key in list(self.p.keys()):
            self.Ps_symbols[key] = symbols(key)
        p=self.Ps_symbols
        if self.setup['rhs']=="RM_forced":
            """ Rosenzweig - MacArthur Model """
            from sympy.functions import sin as symsin
            forcing=(1.0+p['a']*symsin(2.0*np.pi*p['omegaf']*t))
            predation_exp=((p['c']*x*y)/(p['b']+x))
            self.dxdt_eq  = p['r']*x*(1.0-x/p['K'])-predation_exp
            self.dydt_eq  = forcing*p['e']*predation_exp-p['d']*y
        if self.setup['rhs']=="RM_forced_Allee_effect":
            """ Rosenzweig - MacArthur Model with Allee effect """
            from sympy.functions import sin as symsin
            forcing=(1.0+p['a']*symsin(2.0*np.pi*p['omegaf']*t))
            predation_exp=((p['c']*x*y)/(p['b']+x))
            self.dxdt_eq  = p['r']*x*(1.0-x/p['K'])*(x/p['K2']-1.0)-predation_exp
            self.dydt_eq  = forcing*p['e']*predation_exp-p['d']*y
        if self.setup['rhs']=="RM_Allee_effect":
            """ Rosenzweig - MacArthur Model with Allee effect """
            from sympy.functions import sin as symsin
            predation_exp=((p['c']*x*y)/(p['b']+x))
            self.dxdt_eq  = p['r']*x*(1.0-x/p['K'])*(x/p['K2']-1.0)-predation_exp
            self.dydt_eq  = p['e']*predation_exp-p['d']*y
        elif self.setup['rhs']=="RM":
            """ Normal mode """
            from sympy.functions import sin as symsin
            predation_exp=((p['c']*x*y)/(p['b']+x))
            self.dxdt_eq  = p['r']*x*(1.0-x/p['K'])-((p['c']*x*y)/(p['b']+x))
            self.dydt_eq  = p['e']*((p['c']*x*y)/(p['b']+x))-p['d']*y
        """ Creating numpy functions """
        from sympy import solve,Eq
        self.xs_sym=p['b']*p['d']/(p['e']*p['c']-p['d'])
        self.ys_sym=solve(Eq(self.dxdt_eq/x,0),y)[0].subs(x,self.xs_sym)
        self.xs = ufuncify([p['r'],p['K'],p['e'],p['a'],p['omegaf']],[self.sub_parms(self.xs_sym)])
        self.ys = ufuncify([p['r'],p['K'],p['e'],p['a'],p['omegaf']],[self.sub_parms(self.ys_sym)])
        symeqs = Matrix([self.dxdt_eq,self.dydt_eq])
        self.ode  = lambdify((x,y,t,p['r'],p['K'],p['e'],p['a'],p['omegaf']),self.sub_parms(symeqs),"numpy",dummify=False)
        self.dxdt = ufuncify([x,y,t,p['r'],p['K'],p['e'],p['a'],p['omegaf']],[self.sub_parms(self.dxdt_eq)])
        self.dydt = ufuncify([x,y,t,p['r'],p['K'],p['e'],p['a'],p['omegaf']],[self.sub_parms(self.dydt_eq)])
        self.predation_eq = ufuncify([x,y,t,p['r'],p['K'],p['e'],p['a'],p['omegaf']],[self.sub_parms(predation_exp)])
        localJac   = symeqs.jacobian(Matrix([x,y]))
        self.sym_localJac = localJac
        self.localJac = lambdify((x,y,t,p['r'],p['K'],p['e'],p['a'],p['omegaf']),self.sub_parms(localJac),"numpy",dummify=False)
        if self.setup['setPDE'] and self.setup['analyze']:
            self.dbdb  = ufuncify([x,y,p['chi'],p['beta']],[self.sub_parms(localJac[0,0])])
            self.dbds1 = ufuncify([x,y,p['chi'],p['beta']],[self.sub_parms(localJac[0,1])])
            self.dbds2 = ufuncify([x,y,p['chi'],p['beta']],[self.sub_parms(localJac[0,2])])
            self.ds1db = ufuncify([x,y,p['chi'],p['beta']],[self.sub_parms(localJac[1,0])])
            k = symbols('k')
            delta_y  = symbols('delta_y')
            symeqs_lin_analysis = Matrix([self.dxdt_eq-x*k*k,self.dydt_eq-y*delta_y*k*k])
            jaclinanalysis = symeqs_lin_analysis.jacobian(Matrix([x,y]))
            self.symbolic_jaclinanalysis = jaclinanalysis
            self.jaclinanalysis = lambdify((x,y,k),self.sub_parms(jaclinanalysis),"numpy",dummify=False)
        if self.verbose:
            self.print_equations()
            print("Local Jacobian:" ,localJac)
            if self.setup['setPDE'] and self.setup['analyze']:
                print("Linear analysis Jacobian: ", jaclinanalysis)

    """ Printing and parameters related functions """
    def print_parameters(self):
        print(self.p)
    def print_equations(self,numeric=False):
        if numeric:
            print("dxdt = ", self.sub_all_parms(self.dxdt_eq))
            print("dydt = ", self.sub_all_parms(self.dydt_eq))
        else:
            print("dxdt = ", self.dxdt_eq)
            print("dydt = ", self.dydt_eq)
    def print_latex_equations(self):
        from sympy import latex
        print("\partial_t x = ",latex(self.dxdt_eq))
        print("\partial_t y = ",latex(self.dydt_eq))
    """ Functions for use with scipy methods """
    def calc_ode_eigs(self,state,t=0,**kwargs):
        x,y = state[0],state[1],state[2]
        if kwargs:
            self.update_parameters(kwargs)
        return linalg.eigvals(self.localJac(x,y,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf']))
    def calc_Floquet_multipliers(self,state,start=0,step=0.1,method='BDF',**kwargs):
        """ Algorithm taken from https://doi.org/10.1007/s12080-008-0016-2 """
        if kwargs:
            self.update_parameters(kwargs)
        yr = self.p['conv_T_to_t']
        X = np.eye(self.setup['nvar'])
        t = np.arange(start,start+step, step)
        sol_array=solve_ivp(fun=self.scipy_ode_rhs,t_span=(t[0],t[-1]),
                      y0=state,method=method,
                      t_eval=t,jac=self.scipy_ode_jac)
        sol = sol_array.y.T
        time = sol_array.t
        for i,ti in enumerate(time[:-1]):
            state = sol[i,:]
            x,y = state[0],state[1],state[2]
            At = self.localJac(x,y,ti,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'],self.p['omegaf'])
            X = X + step*np.matmul(At,X)
        return linalg.eigvals(X)
    def calc_SpatODE_eigs(self,x,y):
        return linalg.eigvals(self.SpatODEjac(x,y))

    def sigma_k_scan(self,x,y,k_range=[0,1.0],n=1000):
        k_range = np.linspace(k_range[0],k_range[1],n)
        MaxReSigma = np.zeros(n)
        MaxImSigma = np.zeros(n)
        for i,k in enumerate(k_range):
            eigsvalues=linalg.eigvals(self.jaclinanalysis(x,y,k))
            MaxReSigma[i]=np.amax(np.real(eigsvalues))
            MaxImSigma[i]=np.imag(eigsvalues[np.argmax(np.real(eigsvalues))])
        return np.array([k_range,MaxReSigma,MaxImSigma])

    def calc_linear_stability_analysis(self,x,y,k_range=[0,1.0],n=1000):
        k_scan = self.sigma_k_scan(x,y,k_range=[0,0.1],n=1000)
        return k_scan[0][np.argmax(k_scan[1])],np.amax(k_scan[1])
    """ Utilities """
    def sub_parms(self,eqs):
        x,y,t = symbols('x y t')
        for key in list(self.p.keys()):
#            print key
            if key!='r' and key!='K' and key!='e' and key!='a' and key!='omegaf':
                eqs=eqs.subs(self.Ps_symbols[key],self.p[key])
        return eqs
    def sub_all_parms(self,eqs):
        x,y,t = symbols('x y t')
        for key in list(self.p.keys()):
            eqs=eqs.subs(self.Ps_symbols[key],self.p[key])
        return eqs

    """ Spatial functions """
    def set_integrator(self):
        if 'stochastic' not in list(self.setup.keys()):
            integrator_type = {}
            integrator_type['euler'] = self.euler_integrate
            integrator_type['scipy'] = self.ode_integrate
            integrator_type['rk4'] = self.rk4_integrate
            integrator_type['pseudo_spectral'] = self.pseudo_spectral_integrate
            try:
                self.integrator = integrator_type[self.setup['it']]
            except KeyError:
                raise  ValueError("No such integrator : %s."%self.setup['it'])
            if self.setup['it']=='pseudo_spectral':
                self.dt*=100.0
        elif 'stochastic' in list(self.setup.keys()):
#            self.setup_stochastic_rainfall()
            integrator_type = {}
            integrator_type['euler'] = self.euler_integrate_stochastic
            integrator_type['rk4'] = self.rk4_integrate_stochastic
            integrator_type['scipy'] = self.scipy_integrate_stochastic
            try:
                self.integrator = integrator_type[self.setup['it']]
            except KeyError:
                raise  ValueError("No such integrator : %s."%self.setup['it'])

    def rhs_pde(self,state,t=0):
        x,y=np.split(state,3)
        return np.ravel((self.dxdt(x,y,t,self.p['r'],self.p['K'],self.p['e'],self.p['a'],self.p['omegaf']),
                         self.dydt(x,y,t,self.p['r'],self.p['K'],self.p['e'],self.p['a'],self.p['omegaf']))) + self.lapmat*state

    def rhs_ode(self,state,t=0):
        x,y=state
        return self.ode(x,y,t,self.p['r'],self.p['K'],self.p['e'],self.p['a'],self.p['omegaf']).T[0]
    def predation(self,state,t=0):
        x,y=state
        return self.predation_eq(x,y,t,self.p['r'],self.p['K'],self.p['e'],self.p['a'],self.p['omegaf']).T[0]
    def scipy_ode_rhs(self,t,state):
        x,y=state
        return np.squeeze(self.ode(x,y,t,self.p['r'],self.p['K'],self.p['e'],self.p['a'],self.p['omegaf']))
    def scipy_ode_rhs_stochastic(self,t,state):
        print("ME")
        self.add_stochasticity(t)
        x,y=state
        return np.squeeze(self.ode(x,y,t,self.p['r'],self.p['K'],self.p['e'],self.p['a'],self.p['omegaf']))
    def scipy_ode_jac(self,t,state):
        x,y=state
        return self.localJac(x,y,t,self.p['r'],self.p['K'],self.p['e'],self.p['a'],self.p['omegaf'])
    def calc_pde_analytic_jacobian(self,state):
        x,y=np.split(state,3)
        dbdb= sparse.diags(self.dbdb(x,y,self.p['chi'],self.p['beta']))
        dbdw= sparse.diags(self.dbds1(x,y,self.p['chi'],self.p['beta']))
        dbdh= sparse.diags(self.dbds2(x,y,self.p['chi'],self.p['beta']))
        dwdb= sparse.diags(self.ds1db(x,y,self.p['chi'],self.p['beta']))
        dwdw= sparse.diags(self.ds1ds1(x,y,self.p['chi'],self.p['beta']))
        dwdh= sparse.diags(self.ds1ds2(x,y,self.p['chi'],self.p['beta']))
        dhdb= sparse.diags(self.ds2db(x,y,self.p['chi'],self.p['beta']))
        dhdw= sparse.diags(self.ds2ds1(x,y,self.p['chi'],self.p['beta']))
        dhdh= sparse.diags(self.ds2ds2(x,y,self.p['chi'],self.p['beta']))
        local  = sparse.bmat([[dbdb,dbdw,dbdh],
                              [dwdb,dwdw,dwdh],
                              [dhdb,dhdw,dhdh]])
        return sparse.csc_matrix(local)+sparse.csc_matrix(self.lapmat)

    def calc_ode_numerical_jacobian(self,x,y,delta=0.00000001):
        state = np.array([x,y])
        jacobian = []
        for j in range(len(state)):
            state_plus = np.copy(state)
            state_minus = np.copy(state)
            state_plus[j] = state_plus[j]+delta
            state_minus[j] = state_minus[j]-delta
            jacobian.append((np.array(self.dudt(state_plus))-np.array(self.dudt(state_minus)))/(2.0*delta))
        return np.array(jacobian).T
    def check_pde_jacobians(self,n=100):
        import time
        timer_analytic=0
        timer_numeric=0
        error = np.zeros(n)
        for i in range(n):
            print(i)
            x=np.random.random(self.setup['n'])
            y=np.random.random(self.setup['n'])
            state=np.ravel((x,y))
            start_time=time.time()
            numeric=self.calc_pde_numerical_jacobian(state)
            mid_time=time.time()
            analytic=self.calc_pde_analytic_jacobian(state)
            end_time=time.time()
            timer_numeric+=(mid_time-start_time)
            timer_analytic+=(end_time-mid_time)
            error[i]=np.max(np.abs(numeric-analytic))
        print("Max difference is ",np.max(error), ", and mean difference is ",np.mean(error))
        print("Average speed for numeric ", timer_numeric/float(n))
        print("Average speed for analytic ", timer_analytic/float(n))
        print("Analytic ", float(timer_numeric)/float(timer_analytic)," faster.")

    def calc_pde_numerical_jacobian(self,state,delta=0.00000001):
        n = len(state)
        jacobian = []
        for j in range(n):
            state_plus = np.copy(state)
            state_minus = np.copy(state)
            state_plus[j] = state_plus[j]+delta
            state_minus[j] = state_minus[j]-delta
            jacobian.append((self.rhs_pde(state_plus)-self.rhs_pde(state_minus))/(2.0*delta))
        return np.array(jacobian).T

    def calc_numeric_pde_eigs(self,state):
        return linalg.eigvals(self.calc_pde_numerical_jacobian(state))
    def calc_analytic_pde_eigs(self,state):
        return sparse.linalg.eigs(self.calc_pde_analytic_jacobian(state),k=3)[0]

    def check_convergence(self,state,previous_state,tolerance=1.0e-5):
        return np.max(np.abs(state-previous_state))<tolerance
    """Stochastic rainfall functions   """
    def setup_stochastic_rainfall(self,start):
        if self.verbose:
            print("Setting up stochastic term")
#        if kwargs:
#            self.update_parameters(kwargs)
        self.p0 = self.p['p']
        self.years=start+3.0*self.p['conv_T_to_t']/4.0
        self.p_array=[self.p['p']]
        self.years_array=[self.years]
#        self.lim=self.p0*self.dp
        self.dist_name=self.setup['stochastic'][0]
        if type(self.setup['stochastic'][1])==str:
            self.dist_parms=np.loadtxt(self.setup['stochastic'][1])
        else:
            self.dist_parms=self.setup['stochastic'][1]
        distribution_type = {}
        distribution_type['norm_around_p'] = self.norm_around_p
        distribution_type['pearson3_around_p'] = self.pearson3_around_p
        distribution_type['norm'] = self.norm
        distribution_type['chi'] = self.chi
        try:
            self.random_p = distribution_type[self.setup['stochastic'][0]]
        except KeyError:
            raise  ValueError("No such probability distribution : %s."%self.dist_name)
        np.random.seed(self.setup['stochastic'][2])

    def chi(self):
        p=chi.rvs(*self.dist_parms[:-2], loc=self.dist_parms[-2], scale=self.dist_parms[-1])
        return max(0,p)
    def norm(self):
        p=norm.rvs(loc=self.dist_parms[-2], scale=self.dist_parms[-1])
        return max(0,p)
    def norm_around_p(self):
        p=norm.rvs(loc=self.p0, scale=self.dist_parms)
        return max(0,p)
    def pearson3_around_p(self):
        p=pearson3.rvs(0,loc=self.p0, scale=self.dist_parms)
        return max(0,p)

    def create_stochastic_p_time_series(self,time):
        p_time_series = np.zeros_like(time)
        years=time[0]+3.0*self.p['conv_T_to_t']/4.0
        j=0
        for i,t in enumerate(time):
            if t-years>=self.p['conv_T_to_t']:
                years+=self.p['conv_T_to_t']
                j+=1
            p_time_series[i]=self.p_t(t,self.p_array[j],self.p['a'],self.p['omegaf'])
        return p_time_series

    def add_stochasticity(self,time):
        if time-self.years>=self.p['conv_T_to_t']:
            self.years+=self.p['conv_T_to_t']
            self.years_array.append(self.years)
#            lim = self.p0*self.dp
#            X=np.random.normal(scale=self.scale)
#            if np.fabs(X)>lim:
#                X=lim*np.sign(X)
            self.p['p']=self.random_p()
#            self.p['p']=self.p0+X
            self.p_array.append(self.p['p'])
            if self.verbose:
                print("year=",self.years/self.p['conv_T_to_t'],"p=",self.p['p'])
            
    def update_parameters(self,parameters):
        intersection=[i for i in list(self.p.keys()) if i in parameters]
        if intersection:
            for key in intersection:
                if self.setup['verbose']:
                    print(str(key)+"="+str(parameters[key]))
                self.p[key]=parameters[key]

    """Generic integration function                """
    def integrate(self,initial_state=None,step=10,
                  max_time = 1000,tol=1.0e-5,plot=False,savefile=None,
                  create_movie=False,check_convergence=True,
                  sim_step=None,**kwargs):
        if kwargs:
            self.update_parameters(kwargs)
        self.filename = savefile
        if initial_state is None:
            initial_state = self.initial_state
        self.time_elapsed=0
        if 'stochastic' in list(self.setup.keys()):
            self.setup_stochastic_rainfall()
        if sim_step is None:
            self.sim_step=0
        else:
            self.sim_step=sim_step
        if savefile is not None:
            hn.setup_simulation(savefile,self.p,self.setup)
            hn.save_sim_snapshot(savefile,self.sim_step,self.time_elapsed,
                                 self.split_state(initial_state),self.setup)
#        old_result = initial_state
        converged=False
        result = []
        result.append(initial_state)
#        t,result = self.integrator(initial_state,p=p,step=10,finish=10,savefile=self.filename)
        if self.setup['verbose']:
            start=time.time()
            print("Step {:4d}, Time = {:5.1f}".format(self.sim_step,self.time_elapsed))
        while not converged and self.time_elapsed<=max_time:
            old_result = result[-1]
            t,result = self.integrator(result[-1],step=step,finish=step)
#            self.time_elapsed+=t[-1]
            self.sim_step=self.sim_step+1
            if savefile is not None:
                hn.save_sim_snapshot(savefile,self.sim_step,self.time_elapsed,
                                     self.split_state(result[-1]),self.setup)            
            if self.setup['verbose']:
                print("Step {:4d}, Time = {:10.6f}, diff = {:7f}".format(self.sim_step,self.time_elapsed,np.max(np.abs(result[-1]-old_result))))
            if check_convergence:
                converged=self.check_convergence(result[-1],old_result,tol)
                if converged:
                    print("Convergence detected")
        if self.setup['verbose']:
            print("Integration time was {} s".format(time.time()-start))
        if savefile is not None and create_movie:
            print("Creating movie...")
            hn.create_animation_b(savefile)
        return result[-1]

    def calc_stable_state(self,perturb=0.001):
        xs=self.xs(self.p['r'],self.p['K'],self.p['e'],self.p['a'],self.p['omegaf'])
        ys=self.ys(self.p['r'],self.p['K'],self.p['e'],self.p['a'],self.p['omegaf'])
        return [xs-perturb,ys+perturb]
    """ Integrators step functions """
    def ode_integrate(self,initial_state,step=0.01,start=0,finish=100,
                          method='BDF',**kwargs):
        """ Using the new scipy interface to BDF method for stiff equations
        with option to switch to other methods
        """
        if kwargs:
            self.update_parameters(kwargs)
        if type(initial_state)==str:
            if initial_state=='stable_state':
                initial_state=self.calc_stable_state()
        t = np.arange(start,finish+step, step)
        if method=='BDF':
            sjac=self.scipy_ode_jac
        else:
            sjac=None
        sol=solve_ivp(fun=self.scipy_ode_rhs,t_span=(t[0],t[-1]),
                      y0=initial_state,method=method,max_step=step/10.0,
                      t_eval=t,jac=sjac)
        return sol.t,sol.y.T

    def euler_integrate(self,initial_state=None,step=0.1,finish=1000,**kwargs):
        """ Integration using Foreward Euler step
        """
        if kwargs:
            self.update_parameters(kwargs)
        if initial_state is None:
            initial_state=self.state
        time = np.arange(0,finish+step,step)
        result=np.zeros((len(time),len(initial_state)))
        t=0
        result[0]=initial_state
        for i,tout in enumerate(time[1:]):
            old=result[i]
            while t < tout:
                new=old+self.dt*self.rhs(old,self.time_elapsed)
                old=new
                t+=self.dt
                self.time_elapsed+=self.dt
            result[i+1]=old
        self.state=result[-1]
        return time,result
    def euler_integrate_stochastic(self,initial_state=None,step=0.1,finish=1000,**kwargs):
        """ Integration using Foreward Euler step
        """
        if kwargs:
            self.update_parameters(kwargs)
        if initial_state is None:
            initial_state=self.state
        time = np.arange(0,finish+step,step)
        self.p0=self.p['p']
        result=np.zeros((len(time),len(initial_state)))
        t=0
        result[0]=initial_state
        for i,tout in enumerate(time[1:]):
            old=result[i]
            while t < tout:
                self.add_stochasticity(t)
                new=old+self.dt*self.rhs(old,self.time_elapsed)
                old=new
                t+=self.dt
                self.time_elapsed+=self.dt
            result[i+1]=old
        self.state=result[-1]
        self.p_array=np.array(self.p_array)
        return time,result
    
    def rk4_integrate(self,initial_state=None,step=0.1,finish=1000,**kwargs):
        """ Integration using The Runge–Kutta method of 4th order as in:
        https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method
        """
        if kwargs:
            self.update_parameters(kwargs)
        if initial_state is None:
            initial_state=self.state
        time = np.arange(0,finish+step,step)
        result=np.zeros((len(time),len(initial_state)))
        t=0
        result[0]=initial_state
        for i,tout in enumerate(time[1:]):
            old=result[i]
            while t < tout:
                k1=self.rhs(old,self.time_elapsed)
                k2=self.rhs(old+0.5*self.dt*k1,self.time_elapsed+(self.dt/2.0))
                k3=self.rhs(old+0.5*self.dt*k2,self.time_elapsed+(self.dt/2.0))
                k4=self.rhs(old+self.dt*k3,self.time_elapsed+(self.dt))
                new=old+(self.dt/6.0)*(k1+2.0*k2+2.0*k3+k4)
                old=new
                t+=self.dt
                self.time_elapsed+=self.dt
            result[i+1]=old
        self.state=result[-1]
        return time,result

    def scipy_integrate_stochastic(self,initial_state=None,
                                   step=0.1,start=0,finish=1000,
                                   method='RK45',**kwargs):
        """ Using the new scipy interface to RK45 method
        with option to switch to other methods
        """
        if kwargs:
            self.update_parameters(kwargs)
        if initial_state is None:
            initial_state=self.state
        t = np.arange(start,finish+step,step)
        self.setup_stochastic_rainfall(start=start)
        sol=solve_ivp(fun=self.scipy_ode_rhs_stochastic,t_span=(t[0],t[-1]),
                      y0=initial_state,method=method,max_step=self.p['conv_T_to_t'],
                      t_eval=t)
        return sol.t,sol.y.T  
     
    def rk4_integrate_stochastic(self,initial_state=None,
                                 step=0.1,start=0,finish=1000,**kwargs):
        """ Integration using The Runge–Kutta method of 4th order as in:
        https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method
        """
        print("Integration using stochastic rk4 step")
        if kwargs:
            self.update_parameters(kwargs)
        if initial_state is None:
            initial_state=self.state
        time = np.arange(start,finish+step,step)
        self.setup_stochastic_rainfall(start=start)
        result=np.zeros((len(time),len(initial_state)))
        t=start
        result[0]=initial_state
        self.prec=np.zeros(len(time))
        self.prec[0]=self.p0
        for i,tout in enumerate(time[1:]):
            old=result[i]
            while t < tout:
                self.add_stochasticity(t)
                k1=self.rhs(old,self.time_elapsed)
                k2=self.rhs(old+0.5*self.dt*k1,self.time_elapsed+(self.dt/2.0))
                k3=self.rhs(old+0.5*self.dt*k2,self.time_elapsed+(self.dt/2.0))
                k4=self.rhs(old+self.dt*k3,self.time_elapsed+(self.dt))
                new=old+(self.dt/6.0)*(k1+2.0*k2+2.0*k3+k4)
                old=new
                t+=self.dt
                self.time_elapsed+=self.dt
            result[i+1]=old
            self.prec[i+1]=self.p_t(t,self.p['p'],self.p['a'],self.p['omegaf'])
        self.state=result[-1]
        self.p_array=np.array(self.p_array)
        return time,result
    def pseudo_spectral_integrate(self,initial_state=None,step=0.1,finish=1000,**kwargs):
#        print "Integration using pseudo-spectral step"
        if kwargs:
            self.update_parameters(kwargs)
        result=[]
        time = np.arange(0,finish+step,step)
        t=0
        result.append(initial_state)
        for tout in time[1:]:
            self.state=result[-1]
            x,y=self.state.reshape(self.setup['nvar'],*self.setup['n'])
            self.fftb=fftn(x)
            self.ffts1=fftn(y)
            self.ffts2=fftn(s2)
            while t < tout:
                self.fftb = self.multb*(self.fftb + self.dt*fftn(self.dxdt(x,y,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'])))#.real
                self.ffts1 = self.mults1*(self.ffts1 + self.dt*fftn(self.dydt(x,y,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'])))#.real
                self.ffts2 = self.mults2*(self.ffts2 + self.dt*fftn(self.ds2dt(x,y,t,self.p['p'],self.p['chi'],self.p['beta'],self.p['a'])))#.real
                x= ifftn(self.fftb).real
                y= ifftn(self.ffts1).real
                s2= ifftn(self.ffts2).real
                t+=self.dt
                self.time_elapsed+=self.dt
            self.state=np.ravel((x,y))
            self.sim_step+=1
            result.append(self.state)
        return time,result

    def spectral_multiplier(self,dt):
        n=self.setup['n']
        nx=n[0]
        dx=self.dx
        # wave numbers
        k=2.0*np.pi*np.fft.fftfreq(nx,dx)
        if len(n)==1:
            k2=k**2
        if len(n)==2:
            k2=np.outer(k,np.ones(nx))**2
            k2+=np.outer(np.ones(nx),k)**2
        # multiplier
        self.multb = np.exp(-dt*k2)
        self.mults1= np.exp(-dt*self.p['delta_s']*k2)
        self.mults2= np.exp(-dt*self.p['delta_s']*k2)
    """ Auxilary root finding functions """
    def ode_root(self,initial_state,p=None,chi=None,beta=None,a=None,omegaf=None):
        """ """
        if p is not None:
            self.p['p']=p
        if chi is not None:
            self.p['chi']=chi
        if beta is not None:
            self.p['beta']=beta
        if a is not None:
            self.p['a']=a
        if a is not None:
            self.p['omegaf']=omegaf
        sol = root_ode(self.rhs_ode, initial_state)
        return sol.x
    def pde_root(self,initial_state,p=None,chi=None,beta=None,a=None,omegaf=None, fixiter=100,tol=6e-6,smaxiter=1000,integrateiffail=False):
        """ """
        if p is not None:
            self.p['p']=p
        if chi is not None:
            self.p['chi']=chi
        if beta is not None:
            self.p['beta']=beta
        if a is not None:
            self.p['a']=a
        if a is not None:
            self.p['omegaf']=omegaf
        try:
            sol = newton_krylov(self.rhs_pde, initial_state,iter=fixiter, method='lgmres', verbose=int(self.setup['verbose']),f_tol=tol,maxiter=max(fixiter+1,smaxiter))
            converged=True
        except NoConvergence:
            converged=False
            if self.setup['verbose']:
                print("No Convergence flag")
            sol=initial_state
            if integrateiffail:
                print("Integrating instead of root finding")
                sol=self.integrate_till_convergence(initial_state,p)
        return sol,converged

    def setup_initial_condition(self,Vs,**kwargs):
        n = self.setup['n']
        if type(Vs)==str:
            if Vs == "random":
                x = np.random.random(n)*0.5 + 0.1
                y= np.random.random(n)*(0.1) + 0.05
                s2= np.random.random(n)*(0.1) + 0.05
            if Vs == "bare":
                x = np.zeros(n)
                y= np.random.random(n)*(0.1) + 0.05
                s2= np.random.random(n)*(0.1) + 0.05
            elif Vs == "tile":
                fields = kwargs.get('fields', None)
                x,y = np.split(fields,self.setup['nvar'])
                x  = np.tile(x,(self.setup['n'][0],1))
                y = np.tile(y,(self.setup['n'][0],1))
                s2 = np.tile(s2,(self.setup['n'][0],1))                
            elif Vs == "half":
#                import matplotlib.pyplot as plt
                p = kwargs.get('p', None)
#                w = kwargs.get('w', None)
                if p==None:
                    p=self.p['p']
                else:
                    self.p['p']=p
                t,result=self.integrate_ode_bdf([0.01,0.2,0.2],p)
#                plt.plot(t,result[0])
                b_0,s1_0,s2_0 = result.T[-1]
                t,result=self.integrate_ode_bdf([0.9,0.2,0.2],p)
#                plt.plot(t,result[0])
                b_s,s1_s,s2_s = result.T[-1]
                x = np.ones(n)*b_0
                y= np.ones(n)*s1_0
                s2= np.ones(n)*s2_0
                half = int(self.setup['n'][0]/2)
#                width_n = int((float(w)/(2.0*self.setup['l'][0]))*self.setup['n'][0])
                x[:half]=b_s
                y[:half]=s1_s
                s2[:half]=s2_s
#                plt.show()
            elif Vs == "halfrandom":
#                import matplotlib.pyplot as plt
                p = kwargs.get('p', None)
#                w = kwargs.get('w', None)
                if p==None:
                    p=self.p['p']
                else:
                    self.p['p']=p
                t,result=self.integrate_ode_bdf([0.01,0.2,0.2],p)
#                plt.plot(t,result[0])
                b_0,s1_0,s2_0 = result.T[-1]
                t,result=self.integrate_ode_bdf([0.9,0.2,0.2],p)
#                plt.plot(t,result[0])
                b_s,s1_s,s2_s = result.T[-1]
                x = np.ones(n)*b_0
                y= np.ones(n)*s1_0
                s2= np.ones(n)*s2_0
                half = int(self.setup['n'][0]/2)
#                width_n = int((float(w)/(2.0*self.setup['l'][0]))*self.setup['n'][0])
                x[:half]=b_s*np.random.random(n)[:half]
                y[:half]=s1_s*np.random.random(n)[:half]
                s2[:half]=s2_s*np.random.random(n)[:half]
#                plt.show()
            self.initial_state = np.ravel((x,y))
        else:
            self.initial_state = Vs
        self.state = self.initial_state
        if self.setup['it'] == 'pseudo_spectral' and self.setup['setPDE']:
            self.spectral_multiplier(self.dt)
    """ Plot functions """
    def plotLastFrame(self,initial_state=None,p=None,chi=None,beta=None,a=None,savefile=None):
        sol=self.integrate_till_convergence(initial_state,p,chi,beta,a,savefile)
        self.plot(sol)
        return sol
    def split_state(self,state):
        return state.reshape(self.setup['nvar'],*self.setup['n'])
    
    def plot(self,state,fontsize=12,update=False):
        import matplotlib.pylab as plt
        if update:
            plt.ion()
        x,y=state.reshape(self.setup['nvar'],*self.setup['n'])
        if len(self.setup['n'])==1:
            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
    #        x = np.linspace(0,self.setup['l'][0],self.setup['n'][0])
            ax1.set_ylim([-0.1,1.0])
            ax2.set_ylim([0.0,self.p['s_fc']])
            ax3.set_ylim([0.0,self.p['s_fc']])
            ax1.plot(self.X,x)
            ax1.set_xlim([0,self.X[-1]])
            ax2.plot(self.X,y)
            ax3.plot(self.X,s2)
            ax1.set_title(r'$x$', fontsize=fontsize)
            ax2.set_title(r'$s_1$', fontsize=fontsize)
            ax3.set_title(r'$s_2$', fontsize=fontsize)
        elif len(self.setup['n'])==2:
            fig, (ax1, ax2,ax3) = plt.subplots(1, 3, sharex=True, sharey=True)
            fig.subplots_adjust(right=0.8)
#            ax1.imshow(x,cmap=plt.cm.YlGn, vmin = bmin, vmax = bmax)
            ax1.imshow(x,cmap=plt.cm.YlGn)
            ax1.set_adjustable('box-forced')
            ax1.autoscale(False)
            ax1.set_title(r'$x$', fontsize=fontsize)
            ax2.imshow(y,cmap=plt.cm.Blues)
            ax2.set_adjustable('box-forced')
            ax2.autoscale(False)
            ax2.set_title(r'$s_1$', fontsize=fontsize)
#            ax3.imshow(y,cmap=plt.cm.Blues, vmin = smin, vmax = smax)
            ax3.imshow(s2,cmap=plt.cm.Blues)
            ax3.set_adjustable('box-forced')
            ax3.autoscale(False)
            ax3.set_title(r'$s_2$', fontsize=fontsize)
#            im4=ax4.imshow(s2,cmap=plt.cm.Blues, vmin = smin, vmax = smax)
#            cbar_ax2 = fig.add_axes([0.85, 0.54, 0.03, 0.35])
#            fig.colorbar(im1, cax=cbar_ax2)
#            plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    main()