#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 18:15:36 2018

@author: ohm
"""

import netCDF4
import numpy as np

def setup_p_chi_soltype_many_seeds(fname,Ps,Es,p_range,chi_range,seed_range,savefftanalytics):
    """
    Opening an netCDF4 file with 3 groups:
    """
    import time as t
    with netCDF4.Dataset("%s.nc"%fname, 'w', format='NETCDF4') as rootgrp:
        print "Configuring netCDF4 file."
        setup = rootgrp.createGroup('setup')
#        parms = rootgrp.createGroup('Ps')
        setattr(setup, 'nd',int(len(Es['n']))) 
        setattr(setup, 'nx',Es['n'][0]) 
        setattr(setup, 'lx',Es['l'][0])
        if int(len(Es['n']))==2:
            setattr(setup, 'ny',Es['n'][1]) 
            setattr(setup, 'ly',Es['l'][1])
#        for k,v in Ps.items():
#            if k!='dimpar':
#                setattr(parms,k,v)
        rootgrp.description = "Scanning p to chi grid for tlm model with Ps ",Ps
        rootgrp.history = "Created " + t.ctime(t.time())
#        rootgrp.createDimension("x", Es['n'][0])
        rootgrp.createDimension('p', len(p_range))
        rootgrp.createDimension('chi', len(chi_range))
        rootgrp.createDimension('seed', len(seed_range))
#        time = rootgrp.createVariable('time', 'f4', ('time',),zlib=True)
#        time.units = "nondim_year"
#        x = rootgrp.createVariable('x', 'f4', ('x',),zlib=True)
        p   = rootgrp.createVariable('p', 'f4', ('p',),zlib=True)
        chi = rootgrp.createVariable('chi', 'f4', ('chi',),zlib=True)
        seed = rootgrp.createVariable('seed', 'i4', ('seed',),zlib=True)
#        x.units = "nondim_m"
        chi.units = "nondim_trade_off_parameter"
        p.units = "nondim_mmtoyear"
#        x[:] = np.linspace(0,Es['l'][0], Es['n'][0])
        chi[:] = chi_range
        seed[:] = seed_range
        p[:] = p_range
        soltype=rootgrp.createVariable('soltype', 'i4', ('chi','p','seed'),zlib=True)
        shape=soltype[:,:,:].shape
        soltype[:,:,:] = np.zeros(shape,dtype=int)
        if savefftanalytics:
            fft2to1ratio_stochastic=rootgrp.createVariable('fft2to1ratio_stochastic', 'f4', ('chi','p','seed'),zlib=True)
            fft4to1ratio_stochastic=rootgrp.createVariable('fft4to1ratio_stochastic', 'f4', ('chi','p','seed'),zlib=True)
            fft8to1ratio_stochastic=rootgrp.createVariable('fft8to1ratio_stochastic', 'f4', ('chi','p','seed'),zlib=True)
            max_around_fft0_5=rootgrp.createVariable('max_around_fft0_5', 'f4', ('chi','p','seed'),zlib=True)
            max_besides_fft0_5=rootgrp.createVariable('max_besides_fft0_5', 'f4', ('chi','p','seed'),zlib=True)
            fft2to1ratio_stochastic[:,:,:]=np.zeros(shape)
            fft4to1ratio_stochastic[:,:,:]=np.zeros(shape)
            fft8to1ratio_stochastic[:,:,:]=np.zeros(shape)
            max_around_fft0_5[:,:,:]=np.zeros(shape)
            max_besides_fft0_5[:,:,:]=np.zeros(shape)
    print "Output: netCDF file was created: ", fname+".nc"
def setup_p_chi_soltype(fname,Ps,Es,p_range,chi_range):
    """
    Opening an netCDF4 file with 3 groups:
    """
    import time as t
    with netCDF4.Dataset("%s.nc"%fname, 'w', format='NETCDF4') as rootgrp:
        print "Configuring netCDF4 file."
        setup = rootgrp.createGroup('setup')
#        parms = rootgrp.createGroup('Ps')
        setattr(setup, 'nd',int(len(Es['n']))) 
        setattr(setup, 'nx',Es['n'][0]) 
        setattr(setup, 'lx',Es['l'][0])
        if int(len(Es['n']))==2:
            setattr(setup, 'ny',Es['n'][1]) 
            setattr(setup, 'ly',Es['l'][1])
#        for k,v in Ps.items():
#            if k!='dimpar':
#                setattr(parms,k,v)
        rootgrp.description = "Scanning p to chi grid for tlm model with Ps ",Ps
        rootgrp.history = "Created " + t.ctime(t.time())
#        rootgrp.createDimension("x", Es['n'][0])
        rootgrp.createDimension('p', len(p_range))
        rootgrp.createDimension('chi', len(chi_range))
#        time = rootgrp.createVariable('time', 'f4', ('time',),zlib=True)
#        time.units = "nondim_year"
#        x = rootgrp.createVariable('x', 'f4', ('x',),zlib=True)
        p   = rootgrp.createVariable('p', 'f4', ('p',),zlib=True)
        chi = rootgrp.createVariable('chi', 'f4', ('chi',),zlib=True)
#        x.units = "nondim_m"
        chi.units = "nondim_trade_off_parameter"
        p.units = "nondim_mmtoyear"
#        x[:] = np.linspace(0,Es['l'][0], Es['n'][0])
        chi[:] = chi_range
        p[:] = p_range
        rootgrp.createVariable('soltype', 'i4', ('chi','p',),zlib=True)
    print "Output: netCDF file was created: ", fname+".nc"
def setup_p_chi_analyze(fname,Ps,Es,p_range,chi_range,data):
    """
    Opening an netCDF4 file with 3 groups:
    """
    import time as t
    with netCDF4.Dataset("%s.nc"%fname, 'w', format='NETCDF4') as rootgrp:
        print "Configuring netCDF4 file."
        setup = rootgrp.createGroup('setup')
#        parms = rootgrp.createGroup('Ps')
        setattr(setup, 'nd',int(len(Es['n']))) 
        setattr(setup, 'nx',Es['n'][0]) 
        setattr(setup, 'lx',Es['l'][0])
        if int(len(Es['n']))==2:
            setattr(setup, 'ny',Es['n'][1]) 
            setattr(setup, 'ly',Es['l'][1])
#        for k,v in Ps.items():
#            if k!='dimpar':
#                setattr(parms,k,v)
        rootgrp.description = "Scanning p to chi grid for tlm model with Ps ",Ps
        rootgrp.history = "Created " + t.ctime(t.time())
#        rootgrp.createDimension("x", Es['n'][0])
        rootgrp.createDimension('p', len(p_range))
        rootgrp.createDimension('chi', len(chi_range))
        rootgrp.createDimension('fft', len(data['fftBt_stochastic']))
        rootgrp.createDimension('frq', len(data['frq']))
        rootgrp.createDimension('t', len(data['t_stochastic']))
        rootgrp.createDimension('b', len(data['Bt_stochastic']))
#        time = rootgrp.createVariable('time', 'f4', ('time',),zlib=True)
#        time.units = "nondim_year"
#        x = rootgrp.createVariable('x', 'f4', ('x',),zlib=True)
        p   = rootgrp.createVariable('p', 'f4', ('p',),zlib=True)
        chi = rootgrp.createVariable('chi', 'f4', ('chi',),zlib=True)
        frq = rootgrp.createVariable('frq', 'f4', ('frq',),zlib=True)
        t   = rootgrp.createVariable('t', 'f4', ('t',),zlib=True)
#        x.units = "nondim_m"
        chi.units = "nondim_trade_off_parameter"
        p.units = "nondim_mmtoyear"
#        x[:] = np.linspace(0,Es['l'][0], Es['n'][0])
        chi[:] = chi_range
        p[:] = p_range
        frq[:] = data['frq']
        t[:] = data['t_stochastic']/data['conv_T_to_t']
        soltype=rootgrp.createVariable('soltype', 'i4', ('chi','p',),zlib=True)
        average_stochastic_time_series=rootgrp.createVariable('average_stochastic_time_series', 'i4', ('chi','p',),zlib=True)
        fft2to1ratio_stochastic=rootgrp.createVariable('fft2to1ratio_stochastic', 'f4', ('chi','p',),zlib=True)
        fft4to1ratio_stochastic=rootgrp.createVariable('fft4to1ratio_stochastic', 'f4', ('chi','p',),zlib=True)
        fft8to1ratio_stochastic=rootgrp.createVariable('fft8to1ratio_stochastic', 'f4', ('chi','p',),zlib=True)
        max_around_fft0_5=rootgrp.createVariable('max_around_fft0_5', 'f4', ('chi','p',),zlib=True)
        max_besides_fft0_5=rootgrp.createVariable('max_besides_fft0_5', 'f4', ('chi','p',),zlib=True)
        shape = soltype[:,:].shape
        soltype[:,:] = np.zeros(shape,dtype=int)
        average_stochastic_time_series[:,:] = np.zeros(shape,dtype=int)
        print "Shape of soltype",shape
        fft2to1ratio_stochastic[:,:]=np.zeros(shape)
        fft4to1ratio_stochastic[:,:]=np.zeros(shape)
        fft8to1ratio_stochastic[:,:]=np.zeros(shape)
        max_around_fft0_5[:,:]=np.zeros(shape)
        max_besides_fft0_5[:,:]=np.zeros(shape)
        fftBt_stochastic=rootgrp.createVariable('fftBt_stochastic', 'f4', ('chi','p','fft',),zlib=True)
        shape_fft = fftBt_stochastic[:,:,:].shape
        fftBt_stochastic[:,:,:]=np.zeros(shape_fft)
        Bt_stochastic=rootgrp.createVariable('Bt_stochastic', 'f4', ('chi','p','b',),zlib=True)
        shape_Bt = Bt_stochastic[:,:,:].shape
        Bt_stochastic[:,:,:]=np.zeros(shape_Bt)
    print "Output: netCDF file was created: ", fname+".nc"

def setup_Tcollapse_scan(fname,Ps,Es,p_range):
    """
    Opening an netCDF4 file with 3 groups:
    """
    import time as t
    with netCDF4.Dataset("%s.nc"%fname, 'w', format='NETCDF4') as rootgrp:
        print "Configuring netCDF4 file."
        setup = rootgrp.createGroup('setup')
#        parms = rootgrp.createGroup('Ps')
        setattr(setup, 'nd',int(len(Es['n']))) 
        setattr(setup, 'nx',Es['n'][0]) 
        setattr(setup, 'lx',Es['l'][0])
        if int(len(Es['n']))==2:
            setattr(setup, 'ny',Es['n'][1]) 
            setattr(setup, 'ly',Es['l'][1])
#        for k,v in Ps.items():
#            if k!='dimpar':
#                setattr(parms,k,v)
        rootgrp.description = "Scanning p to chi grid for tlm model with Ps ",Ps
        rootgrp.history = "Created " + t.ctime(t.time())
        rootgrp.createDimension('p', len(p_range))
        rootgrp.createDimension('time', len(p_range))
        p   = rootgrp.createVariable('p', 'f4', ('p',),zlib=True)
        p.units = "nondim_mmtoyear"
#        x[:] = np.linspace(0,Es['l'][0], Es['n'][0])
        p[:] = p_range
        time=rootgrp.createVariable('Tcollapse', 'f4', ('time',),zlib=True)
        time.units = "yr"
    print "Output: netCDF file was created: ", fname+".nc"

def setup_locate_transition_along_chi(fname,Ps,Es,chi_range,soltype_pmin,soltype_pmax):
    """
    Opening an netCDF4 file with 3 groups:
    """
    import time as t
    with netCDF4.Dataset("%s.nc"%fname, 'w', format='NETCDF4') as rootgrp:
        print "Configuring netCDF4 file."
        setup = rootgrp.createGroup('setup')
#        parms = rootgrp.createGroup('Ps')
        setattr(setup, 'nd',int(len(Es['n']))) 
        setattr(setup, 'nx',Es['n'][0]) 
        setattr(setup, 'lx',Es['l'][0])
        if int(len(Es['n']))==2:
            setattr(setup, 'ny',Es['n'][1]) 
            setattr(setup, 'ly',Es['l'][1])
#        for k,v in Ps.items():
#            if k!='dimpar':
#                setattr(parms,k,v)
        rootgrp.description = "Scanning to locate transition from soltype {} to soltype{} along chi for tlm model with Ps ".format(soltype_pmin,soltype_pmax),Ps
        rootgrp.history = "Created " + t.ctime(t.time())
#        rootgrp.createDimension("x", Es['n'][0])
        rootgrp.createDimension('p', len(chi_range))
        rootgrp.createDimension('chi', len(chi_range))
#        time = rootgrp.createVariable('time', 'f4', ('time',),zlib=True)
#        time.units = "nondim_year"
#        x = rootgrp.createVariable('x', 'f4', ('x',),zlib=True)
        p   = rootgrp.createVariable('p', 'f4', ('p',),zlib=True)
        chi = rootgrp.createVariable('chi', 'f4', ('chi',),zlib=True)
#        x.units = "nondim_m"
        chi.units = "nondim_trade_off_parameter"
        p.units = "nondim_mmtoyear"
#        x[:] = np.linspace(0,Es['l'][0], Es['n'][0])
        chi[:] = chi_range
    print "Output: netCDF file was created: ", fname+".nc"

def save_p_chi_soltype(fname,p_i,chi_i,p,chi,data,savefftanalytics):
    """ 
    """
    with netCDF4.Dataset("%s.nc"%fname, 'a') as rootgrp:
        rootgrp['soltype'][chi_i,p_i] = data['soltype']
        if savefftanalytics:
            rootgrp['fft2to1ratio_stochastic'][chi_i,p_i] = data['fft2to1ratio_stochastic']
            rootgrp['fft4to1ratio_stochastic'][chi_i,p_i] = data['fft4to1ratio_stochastic']
            rootgrp['fft8to1ratio_stochastic'][chi_i,p_i] = data['fft8to1ratio_stochastic']
            rootgrp['max_around_fft0_5'][chi_i,p_i] = data['max_around_fft0_5']
            rootgrp['max_besides_fft0_5'][chi_i,p_i] = data['max_besides_fft0_5']
def save_p_chi_soltype_per_seed(fname,p_i,chi_i,seed_i,data,savefftanalytics):
    """ 
    """
    with netCDF4.Dataset("%s.nc"%fname, 'a') as rootgrp:
        rootgrp['soltype'][chi_i,p_i,seed_i] = data['soltype']
        if savefftanalytics:
            rootgrp['fft2to1ratio_stochastic'][chi_i,p_i,seed_i] = data['fft2to1ratio_stochastic']
            rootgrp['fft4to1ratio_stochastic'][chi_i,p_i,seed_i] = data['fft4to1ratio_stochastic']
            rootgrp['fft8to1ratio_stochastic'][chi_i,p_i,seed_i] = data['fft8to1ratio_stochastic']
            rootgrp['max_around_fft0_5'][chi_i,p_i,seed_i] = data['max_around_fft0_5']
            rootgrp['max_besides_fft0_5'][chi_i,p_i,seed_i] = data['max_besides_fft0_5']
def save_p_chi_analyze(fname,p_i,chi_i,p,chi,data):
    """ 
    """
    with netCDF4.Dataset("%s.nc"%fname, 'a') as rootgrp:
        rootgrp['soltype'][chi_i,p_i] = data['soltype']
        if data['soltype'] > 1:
            rootgrp['fft2to1ratio_stochastic'][chi_i,p_i] = data['fft2to1ratio_stochastic']
            rootgrp['fft4to1ratio_stochastic'][chi_i,p_i] = data['fft4to1ratio_stochastic']
            rootgrp['fft8to1ratio_stochastic'][chi_i,p_i] = data['fft8to1ratio_stochastic']
            rootgrp['max_around_fft0_5'][chi_i,p_i] = data['max_around_fft0_5']
            rootgrp['max_besides_fft0_5'][chi_i,p_i] = data['max_besides_fft0_5']
            rootgrp['fftBt_stochastic'][chi_i,p_i,:] = data['fftBt_stochastic']

def save_Tcollapse(fname,i,Tcollapse):
    """ 
    """
    with netCDF4.Dataset("%s.nc"%fname, 'a') as rootgrp:
        rootgrp['Tcollapse'][i] = Tcollapse

def save_transition(fname,chi,p):
    """ 
    """
    with netCDF4.Dataset("%s.nc"%fname, 'a') as rootgrp:
        chi_range=rootgrp['chi'][:]
        idx=np.searchsorted(chi_range,chi)
        rootgrp['p'][idx] = p

def load_p_chi_info(fname,printroot=False):
    if not fname.endswith('.nc'):
        fname+='.nc'
    with netCDF4.Dataset(fname, 'r') as rootgrp:
        if printroot:
            print rootgrp
        p = rootgrp['p'][:]
        chi = rootgrp['chi'][:]
    return p,chi
def load_p_chi_data(fname,data):
    if not fname.endswith('.nc'):
        fname+='.nc'
    with netCDF4.Dataset(fname, 'r') as rootgrp:
        data = rootgrp[data][:,:]
    return data

def load_p_chi_soltype(fname,manyseeds=False):
    """ Load results from p to a grid scan
    Returns soltype,p,chi
    """
    if not fname.endswith('.nc'):
        fname+='.nc'
    with netCDF4.Dataset(fname, 'r') as rootgrp:
        p = rootgrp['p'][:]
        chi = rootgrp['chi'][:]
        if manyseeds:
            soltype = rootgrp['soltype'][:,:,:]
        else:
            soltype = rootgrp['soltype'][:,:]
    return soltype,p,chi

def conv_nc_to_dict(fname):
    if not fname.endswith('.nc'):
        fname+='.nc'
    data = {}
    with netCDF4.Dataset(fname, 'r') as rootgrp:
        for key in rootgrp.variables.keys():
            if len(rootgrp[key].shape)==3:
                data[key]=rootgrp[key][:,:,:]
            elif len(rootgrp[key].shape)==2:
                data[key]=rootgrp[key][:,:]
            else:
                data[key]=rootgrp[key][:]
    return data
def save_hdf_from_nc(fname):
    from deepdish.io import save
    data=conv_nc_to_dict(fname)
    if fname.endswith('.nc'):
        fname=fname[:-3]
    save(fname+'.hdf5',data)

def plot_p_achi_grid(fname,pstep,astep,t=None):
    import matplotlib.pylab as plt
    b1,b2,w,p,a=load_p_chi_grid(fname,pstep,astep)
    fig,ax = plt.subplots(1,2)
    if t is None:
        t = [0,1]
    ax[0].imshow(b1,extent=[0,1,t[-1],t[0]],aspect='auto',cmap=plt.cm.YlGn,vmin=0,vmax=1.0)
    ax[1].imshow(b2,extent=[0,1,t[-1],t[0]],aspect='auto',cmap=plt.cm.YlGn,vmin=0,vmax=1.0)
    plt.show()