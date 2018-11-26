#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 19:06:04 2017

@author: ohm
"""
import numpy as np
from tlmModel import tlmModel
import utilities.handle_netcdf_p_chi_grid as hn
from scipy.signal import argrelextrema
import time
import deepdish.io as dd
import matplotlib.pyplot as plt

np.seterr(all='ignore')
yatir = "rain/chi_distribution_mean_annual_rainfall_yatir_1964_2016.dat"
Ps='auto/tlm_set25.hdf5'
Es={'rhs':"LambdaMsGamma",
    'n':(1024,),'l':(256.0,),
    'bc':"neumann",'stochastic':["pearson3_around_p",0.1,100],
    'it':"rk4",'dt':0.001,
    'analyze':False,'verbose':False,'setPDE':False}

def plotBtStochastic(data,plotforced=False,size=6,savefile=None):
    if type(data)==str:
        data=dd.load(data)
    fig,ax=plt.subplots(2,2,figsize=(2*4*1.618,2*4))
    if plotforced:
#        ax[1].plot(data['t_forced']/data['conv_T_to_t'],data['p_time_series'],'g')
        ax[0,0].plot(data['t_forced']/data['conv_T_to_t'],data['Bt_forced'],'b')
    ax[1,0].plot(data['t_stochastic']/data['conv_T_to_t'],data['p_time_series'],'g')
    ax[0,0].plot(data['t_stochastic']/data['conv_T_to_t'],data['Bt_stochastic'],'r')
#    ax[0,0]
    if 'T_consecutive_drop' in data.keys():
        ax[0,0].axvline(data['T_consecutive_drop'],color='m',lw=2)
    ax[0,0].set_ylabel(r'$b$')
    ax[1,0].set_ylabel(r'$p$')
    ax[1,0].set_xlabel(r'$T\,[yr]$')
    ax[0,0].set_xlabel(r'$T\,[yr]$')
#    plt.figure()
    ax[0,1].plot(data['frq'][1:],data['fftBt_stochastic'][1:],'r')
    ax[0,1].plot(data['frq'][1:],data['fftBt_forced'][1:],'b')
    ax[0,1].set_xlabel(r'$freq\,[1/yr]$')
    ax[1,1].set_xlabel(r'$p$')
    ax[0,1].set_xlim([0,2])
#    plt.figure()
    dump=ax[1,1].hist(data['p_array'],normed=1)
    plt.tight_layout()
    if savefile is not None:
        p=data['p']
        chi=data['chi']
        dist_seed=data['dist_seed']
#        a=data['a']
        dist_type=data['dist_type']
        dist_scale=data['dist_scale']
        fname = savefile+"_{}_seed{}_scale{:3.2}_p{:3.2f}_chi{:2.1f}".format(dist_type,str(dist_seed),dist_scale,p,chi).replace(".","_")
        plt.savefig(fname+".png")
        plt.savefig(fname+".pdf")
    plt.show()
    

def check_death(Bt,period,threshold):
    """ Return True if the last Bt[-period:] is bellow the threshold
    Bt - time series
    period - the index of the period I would like to check
    threshold - threshold of declaration of biomass death
    """
    if np.all(np.amax(Bt[-int(period):])<threshold):
        return True
    else:
        return False

def find_consecutive_drop_in_biomass(Bt,number_years=3):
    number_of_negative=0
    number_years+=1
    maximas_args = argrelextrema(Bt,np.greater)[0]
    maximas=np.append(np.diff(Bt[maximas_args]),0)
    for i,indx in enumerate(maximas_args):
        if maximas[i]<0:
            number_of_negative+=1
        else:
            number_of_negative=0
        if number_of_negative==number_years:
            break
    return indx

def find_index_of_death(Bt,check_death_period,threshold=1e-3):
    """ Finding the index from thereon the biomass is considered to be dead
    Returning index of Bt array such that Bt[i:i+period] is bellow threshold
    I use numpy roll to prevent going beyond the array limits
    """
    Bsol = Bt.copy()
    period=check_death_period
    i=0
    dead=np.all(np.amax(Bsol[:int(period)])<threshold)
    while not dead:
        i+=1
        Bsol=np.roll(Bt,-i)
        dead=np.all(np.amax(Bsol[:int(period)])<threshold)
        if i==len(Bt):
            break
    return i
def ticks_to_years(ticks,one_year_ticks):
    return float(ticks)/float(one_year_ticks)
def calc_collapse_time(data,threshold=1e-3):
    ticks=find_index_of_death(data['Bt_stochastic'],data['check_death_period'],threshold=threshold)
    return ticks_to_years(ticks,data['one_year_ticks'])


def bisection(pmin,pmax,soltype_pmin,soltype_pmax,chi,a=1.0,tol=0.005,timeit=False):
    if timeit:
        start=time.time()
    p=(pmin+pmax)/2.0
    while (pmax-pmin)/2.0>tol:
        print "Checking for p=",p
        soltype=check_soltype(p,chi,a)
        print "For p=",p," soltype=",soltype
        if soltype==soltype_pmin:
            pmin = p
        elif soltype==soltype_pmax:
            pmax = p
        else:
            raise UserWarning("The bisection did not find soltype matching the ones in the limits of the section.")
        p=(pmin+pmax)/2.0
    if timeit:
        print "Bisection took ",time.time()-start," seconds"
    return p

def locate_transition_along_chi(chi_range,
                                pmin,pmax,
                                soltype_pmin,soltype_pmax,
                                a=1.0,tol=0.005,
                                fname=None,locker=None,verbose=False):
    p=np.zeros_like(chi_range)
    for i,chi in enumerate(chi_range):
        p[i]=bisection(pmin,pmax,soltype_pmin,soltype_pmax,chi)
    if fname is not None:
        if locker is not None:
            locker.acquire()
        for i,chi in enumerate(chi_range):
            hn.save_transition(fname,chi,p[i])
        if locker is not None:
            locker.release()

def apply_async_locate_transition_along_chi_nc(pmin,pmax,a,
                                               soltype_pmin,soltype_pmax,
                                               fname,
                                               Nchi=100,sections=10,
                                               numproc=10,tol=0.005,
                                               verbose=False,send_email=None):
    import multiprocessing as mp
    import time
    if send_email is not None:
        import getpass
        pswd = getpass.getpass('Password:')
    if verbose:
        print "Scanning p_min={},p_max={},a={}".format(pmin,pmax,a)
    start = time.time()
    chi_range = np.linspace(0,1,Nchi)
    split_chi_range = np.split(chi_range,sections)
    hn.setup_locate_transition_along_chi(fname,Ps,Es,chi_range,soltype_pmin,soltype_pmax)
    from utilities.multiprocess import available_cpu_count
    availble_cpus = int(available_cpu_count() - 2)
    numproc=min(numproc,availble_cpus)
    if verbose:
        print "Using",numproc," processors"
    pool = mp.Pool(processes=numproc)
    m = mp.Manager()
    locker = m.Lock()
    for i,chi_section in enumerate(split_chi_range):
        pool.apply_async(locate_transition_along_chi,
                         args = (chi_section,pmin,pmax,soltype_pmin,soltype_pmax,a,tol,fname,locker,verbose))
    pool.close()
    pool.join()
    print "Took ",time.time()-start
    if verbose:
        print "Saved to:",fname
    if send_email is not None:
        try:
            import smtplib
            from socket import gaierror
            server = smtplib.SMTP('smtp.gmail.com', 587)
            server.ehlo()
            server.starttls()
            server.login(send_email, pswd)
            msg = "\r\n".join([
                    "From: {}".format(send_email),
                    "To: {}".format(send_email),
                    "Subject: Apply async calculations finished",
                    "",
                    "Apply async calculations finished and saved to:", fname
                    ])
            server.sendmail(send_email, send_email, msg)
            server.quit()
        except gaierror:
            pass

def analyze_time_series(p,chi,dist_seed=100,dist_scale=0.1,a=1.0,
                        dist_type="norm_around_p",
                        initial_conditions = [0.9,0.2,0.2],
                        death_threshold=1.0e-3,years_to_declare_death=5,
                        T=100.0,centuries=(10,2),
                        Es=Es,Ps=Ps,
                        savefile=None,timer=False,plot=True):
    if timer:
        import time
        start_time=time.time()
    soltype=None
    data={}
    Es['stochastic']=[dist_type,dist_scale,dist_seed]
    m = tlmModel(Es=Es,Ps=Ps,Vs=None)
    yr = m.p['conv_T_to_t']
    Fs = 2.0*T  # sampling rate
    Ts = 1.0/Fs # sampling interval
    Time=np.arange(0,T,Ts)
    n = len(Time) # length of the signal
    k = np.arange(n)
    TT = n/Fs
    frq = k/(TT) # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    idx1yr = (np.abs(frq - 1.0)).argmin()
    delta_frq = frq[idx1yr]-frq[idx1yr-1]
    indx_10delta_frq = max(1,int(0.05/delta_frq))
    idx0_5yr = (np.abs(frq - 0.5)).argmin()
    idx0_25yr = (np.abs(frq - 0.25)).argmin()
    idx0_125yr = (np.abs(frq - 0.125)).argmin()
    # First time integration for constant precipitation
    # to start the periodic time integration from a steady state such to
    # reduce transients, and also to check whether the state is vegetated
    t_const,sol_const=m.ode_integrate(np.array(initial_conditions),
                          p=p,chi=chi,a=0,
                          finish=T*yr)
    bsol_constant=sol_const[:,0]
    data['t_const']=t_const
    data['Bt_constant']=bsol_constant
    data['conv_T_to_t']=yr
    eigenvalues = m.calc_ode_eigs(sol_const[-1])
    ReSigma=np.amax(np.real(eigenvalues))
    ImSigma=np.imag(eigenvalues[np.argmax(np.real(eigenvalues))])
    fftBt_constant = np.absolute((np.fft.fft(bsol_constant)/n)[range(n/2)])
    data['fftBt_constant']=fftBt_constant
    data['frq']=frq
    # Define parameters for checking if biomass is dead
    one_year_ticks = int(len(bsol_constant)/100.0)
    check_death_period = int(one_year_ticks*years_to_declare_death)
    if check_death(bsol_constant,check_death_period,death_threshold):
#        print "HE IS CONSTANTLY DEAD!"
        soltype = 0
        data['collapse_time']=0
    else:
        # Time integration for periodic precipitation
        t_forced,sol_periodic=m.ode_integrate(sol_const[-1],p=p,chi=chi,a=a,start=T*yr,
                              finish=(1.0+centuries[0])*T*yr)
        bsol_forced=sol_periodic[:,0]
        if check_death(bsol_forced,check_death_period,death_threshold):
#            print "HE IS PERIODICALLY DEAD!"
            soltype = 1
            data['collapse_time']=0
        else:
            means=np.mean(bsol_forced[len(bsol_forced)/2:])
            amplitudes=(np.amax(bsol_forced[len(bsol_forced)/2:])-np.amin(bsol_forced[len(bsol_forced)/2:]))/2.0
            maxs=np.amax(bsol_forced[len(bsol_forced)/2:])
            mins=np.amin(bsol_forced[len(bsol_forced)/2:])
            data.update({'A':amplitudes,'mean':means,'max':maxs,'min':mins,
                         'p':p,'chi':chi,'a':a,
                         'ReSigma':ReSigma,'ImSigma':ImSigma})
            fftBt_forced = np.absolute((np.fft.fft(bsol_forced[-len(bsol_constant):])/n)[range(n/2)])
            power_fftBt_1yr = np.trapz(np.fabs(fftBt_forced[idx1yr-indx_10delta_frq:idx1yr+indx_10delta_frq]))
            power_fftBt_0_5yr = np.trapz(np.fabs(fftBt_forced[idx0_5yr-indx_10delta_frq:idx0_5yr+indx_10delta_frq]))
            power_fftBt_0_25yr = np.trapz(np.fabs(fftBt_forced[idx0_25yr-indx_10delta_frq:idx0_25yr+indx_10delta_frq]))
            power_fftBt_0_125yr = np.trapz(np.fabs(fftBt_forced[idx0_125yr-indx_10delta_frq:idx0_125yr+indx_10delta_frq]))
            fft2to1ratio=power_fftBt_0_5yr/power_fftBt_1yr
            fft4to1ratio=power_fftBt_0_25yr/power_fftBt_1yr
            fft8to1ratio=power_fftBt_0_125yr/power_fftBt_1yr
            data['fft2to1ratio_forced']=fft2to1ratio
            data['fft4to1ratio_forced']=fft4to1ratio
            data['fft8to1ratio_forced']=fft8to1ratio
            data['indx_10delta_frq']=indx_10delta_frq
            data['delta_frq']=delta_frq
            data['idx1yr']=idx1yr
            data['idx0_5yr']=idx0_5yr
            data['idx0_25yr']=idx0_25yr
            data['idx0_125yr']=idx0_125yr
            data['fftBt_forced']=fftBt_forced
            data['Bt_forced']=bsol_forced
            data['t_forced']=t_forced
            # time integration for periodic precipitation with stochastic
            # mean annual precipitation
            t_stochastic,sol_stochastic=m.scipy_integrate_stochastic(sol_periodic[-1],p=p,chi=chi,a=a,
                                                          start=(1.0+centuries[0])*T*yr,
                                                          finish=(1.0+centuries[0]+centuries[1])*T*yr)
            bsol_stochastic=sol_stochastic[:,0]
            if check_death(bsol_stochastic,check_death_period,death_threshold):
#                print "HE IS STOCHASTICALLY DEAD!"
                soltype=2
                indx_consecutive_drop = find_consecutive_drop_in_biomass(bsol_stochastic)
                data['T_consecutive_drop']=(1.0+centuries[0])*T+ticks_to_years(indx_consecutive_drop,one_year_ticks)
            else:
#                print "HE IS STOCHASTICALLY ALIVE!"
                soltype=3
                data['T_consecutive_drop']=(1.0+centuries[0]+centuries[1])*T
            data['Bt_stochastic']=bsol_stochastic
            data['t_stochastic']=t_stochastic
            data['p_array']=np.array(m.p_array)
            data['p_time_series']=m.create_stochastic_p_time_series(t_stochastic)
            fftBt = np.absolute((np.fft.fft(bsol_stochastic[:len(bsol_constant)])/n)[range(n/2)])
            power_fftBt_1yr = np.trapz(np.fabs(fftBt[idx1yr-indx_10delta_frq:idx1yr+indx_10delta_frq]))
            power_fftBt_0_5yr = np.trapz(np.fabs(fftBt[idx0_5yr-indx_10delta_frq:idx0_5yr+indx_10delta_frq]))
            power_fftBt_0_25yr = np.trapz(np.fabs(fftBt[idx0_25yr-indx_10delta_frq:idx0_25yr+indx_10delta_frq]))
            power_fftBt_0_125yr = np.trapz(np.fabs(fftBt[idx0_125yr-indx_10delta_frq:idx0_125yr+indx_10delta_frq]))
            fft2to1ratio=power_fftBt_0_5yr/power_fftBt_1yr
            fft4to1ratio=power_fftBt_0_25yr/power_fftBt_1yr
            fft8to1ratio=power_fftBt_0_125yr/power_fftBt_1yr
            data['fft2to1ratio_stochastic']=fft2to1ratio
            data['fft4to1ratio_stochastic']=fft4to1ratio
            data['fft8to1ratio_stochastic']=fft8to1ratio
            data['max_around_fft0_5'] = np.amax(fftBt[(data['idx0_5yr']-10*data['indx_10delta_frq']):(data['idx0_5yr']+10*data['indx_10delta_frq'])])
            first_half_max = np.amax(fftBt[1:(data['idx0_5yr']-10*data['indx_10delta_frq'])])
            second_half_max = np.amax(fftBt[(data['idx0_5yr']+10*data['indx_10delta_frq']):data['idx1yr']])
            data['max_besides_fft0_5'] = np.amax([first_half_max,second_half_max])
            data['fftBt_stochastic']=fftBt
            data['check_death_period']=check_death_period
            data['death_threshold']=death_threshold
            data['one_year_ticks']=one_year_ticks
    data['soltype']=soltype
    data['p']=p
    data['chi']=chi
    data['dist_seed']=dist_seed=100
    data['a']=a
    data['dist_type']=dist_type
    data['dist_scale']=dist_scale
    if timer:
        print "Analyze took ",str(time.time()-start_time)
    if plot:
        plotBtStochastic(data,savefile=savefile)
    if savefile is not None:
        fname = savefile+"_{}_seed{}_scale{:3.2}_p{:3.2f}_chi{:2.1f}".format(dist_type,str(dist_seed),dist_scale,p,chi).replace(".","_")
        dd.save(fname+".hdf5",data)
    return data

def save_analyze(p,p_i,chi,chi_i,a,dist_type,dist_scale,
                 Es,Ps,
                 average_stochastic_time_series,
                 fname,locker=None,verbose=False):
    if verbose:
        print "Starting p={},chi={}".format(p,chi)
    data = analyze_time_series(p,chi,a,dist_type,dist_scale,
                               average_stochastic_time_series=average_stochastic_time_series,
                               Es=Es,Ps=Ps)
    if locker is not None:
        locker.acquire()
    hn.save_p_chi_analyze(fname,p_i,chi_i,p,chi,data)
    if locker is not None:
        locker.release()

def apply_async_tlm_p_chi_scan_analyze(pmin,pmax,a,
                                       dist_type,dist_scale,
                                       fname,
                                       Nchi=11,dp=0.01,Es=Es,Ps=Ps,
                                       average_stochastic_time_series=0,
                                       numproc=10,savehdf=False,
                                       verbose=False,send_email=None,above_diagonal=None):
    import multiprocessing as mp
    import time
    if send_email is not None:
        import getpass
        pswd = getpass.getpass('Password:')
    start = time.time()
    p_range = np.arange(pmin,pmax+dp,dp)
    chi_range = np.linspace(0,1,Nchi)
    if above_diagonal is not None and len(above_diagonal)==2:
        diagb = float(above_diagonal[1])
        diaga = float(above_diagonal[0])-float(above_diagonal[1])
    if verbose:
        print "Running one calculation to setup the netCDF4 file"
    data = analyze_time_series(pmax,0,a,dist_type,dist_scale,Es=Es)
    if verbose:
        print "Setting up netCDF4 file"
    hn.setup_p_chi_analyze(fname,Ps,Es,p_range,chi_range,data)
    from utilities.multiprocess import available_cpu_count
    availble_cpus = int(available_cpu_count() - 2)
    numproc=min(numproc,availble_cpus)
    if verbose:
        print "Using",numproc," processors"
    pool = mp.Pool(processes=numproc)
    m = mp.Manager()
    locker = m.Lock()
    counter = 0 
    if verbose:
        print "Scanning p_min={},p_max={},a={}".format(pmin,pmax,a)
    for chi_i,chi in enumerate(chi_range):
        for p_i,p in enumerate(p_range):
            if above_diagonal is None:
                counter+=1
                pool.apply_async(save_analyze,
                                 args = (p,p_i,chi,chi_i,a,dist_type,dist_scale,
                                         Es,Ps,
                                         average_stochastic_time_series,
                                         fname,locker,verbose))
            elif above_diagonal is not None and len(above_diagonal)==2:
                pdiag = diaga*chi+diagb
                if p>=pdiag:
                    counter+=1
                    pool.apply_async(save_analyze,
                                     args = (p,p_i,chi,chi_i,a,dist_type,dist_scale,
                                             Es,Ps,
                                             average_stochastic_time_series,
                                             fname,locker,verbose))
            else:
                pass
    pool.close()
    pool.join()
    calc_time = time.time()-start
    if verbose:
        print "Took ",calc_time
        print "Made ", counter," calculations"
        print "Saved to:",fname
    if savehdf:
        from utilities.handle_netcdf_p_chi_grid import load_p_chi_info,load_p_chi_data
        p,chi=load_p_chi_info(fname)
        data = {}
        data['p']=p
        data['chi']=chi
        data['soltype']=load_p_chi_data(fname,"soltype")
        data['averaged']=load_p_chi_data(fname,"average_stochastic_time_series")
        data['max_around_fft0_5']=load_p_chi_data(fname,"max_around_fft0_5")
        data['max_besides_fft0_5']=load_p_chi_data(fname,"max_besides_fft0_5")
        data['fft2to1ratio_stochastic']=load_p_chi_data(fname,"fft2to1ratio_stochastic")
        data['fft4to1ratio_stochastic']=load_p_chi_data(fname,"fft4to1ratio_stochastic")
        data['fft8to1ratio_stochastic']=load_p_chi_data(fname,"fft8to1ratio_stochastic")
        dd.save(fname+".hdf5",data)
    if send_email is not None:
        try:
            import smtplib
            from socket import gaierror
            server = smtplib.SMTP('smtp.gmail.com', 587)
            server.ehlo()
            server.starttls()
            server.login(send_email, pswd)
            msg = "\r\n".join([
                    "From: {}".format(send_email),
                    "To: {}".format(send_email),
                    "Subject: Apply async of analyze to tlm model calculations finished",
                    "",
                    "Apply async calculations finished in ",str(calc_time)," sec, made ", counter," calculations and saved to:", fname
                    ])
            server.sendmail(send_email, send_email, msg)
            server.quit()
        except gaierror:
            pass

def check_soltype(p,chi,a,
                  initial_conditions = [0.9,0.2,0.2],
                  death_threshold=1.0e-3,years_to_declare_death=5,
                  T = 100,savefftanalytics=False,
                  Ps=Ps,Es=Es):
    soltype=None
    m = tlmModel(Es=Es,Ps=Ps,Vs=None)
    yr = m.p['conv_T_to_t']
    data={}
    # First time integration for constant precipitation
    # to start the periodic time integration from a steady state such to
    # reduce transients, and also to check whether the state is vegetated
    t_constant,sol=m.ode_integrate(np.array(initial_conditions),
                          p=p,chi=chi,a=0,
                          finish=T*yr)
    bsol_constant=sol[:,0]
    if savefftanalytics:
        Fs = 2.0*T  # sampling rate
        Ts = 1.0/Fs # sampling interval
        Time=np.arange(0,T,Ts)
        n = len(Time) # length of the signal
        k = np.arange(n)
        TT = n/Fs
        frq = k/(TT) # two sides frequency range
        frq = frq[range(n/2)] # one side frequency range
        idx1yr = (np.abs(frq - 1.0)).argmin()
        delta_frq = frq[idx1yr]-frq[idx1yr-1]
        indx_10delta_frq = max(1,int(0.05/delta_frq))
        idx0_5yr = (np.abs(frq - 0.5)).argmin()
        idx0_25yr = (np.abs(frq - 0.25)).argmin()
        idx0_125yr = (np.abs(frq - 0.125)).argmin()   
        data['indx_10delta_frq']=indx_10delta_frq
        data['idx1yr']=idx1yr
        data['idx0_5yr']=idx0_5yr
        data['idx0_25yr']=idx0_25yr
        data['idx0_125yr']=idx0_125yr
        data['fft2to1ratio_stochastic']=0.0
        data['fft4to1ratio_stochastic']=0.0
        data['fft8to1ratio_stochastic']=0.0
        data['max_around_fft0_5']=0.0
        data['max_besides_fft0_5']=0.0
    one_year_ticks = int(len(bsol_constant)/100.0)
    check_death_period = int(one_year_ticks*years_to_declare_death)
    if check_death(bsol_constant,check_death_period,death_threshold):
#        print "HE IS CONSTANTLY DEAD!"
        soltype = 0
    else:
        # Time integration for periodic precipitation
        t_forced,sol=m.ode_integrate(sol[-1],p=p,chi=chi,a=a,start=T*yr,
                              finish=10*T*yr)
        bsol_forced=sol[:,0]
        if check_death(bsol_forced,check_death_period,death_threshold):
#            print "HE IS PERIODICALLY DEAD!"
            soltype = 1
        else:
            # time integration for periodic precipitation with stochastic
            # mean annual precipitation
            t,sol=m.scipy_integrate_stochastic(sol[-1],p=p,chi=chi,a=a,
                                               start=10*T*yr,
                                               finish=11*T*yr)
            bsol_stochastic=sol[:,0]
            if check_death(bsol_stochastic,check_death_period,death_threshold):
#                print "HE IS STOCHASTICALLY DEAD!"
                soltype=2
            else:
#                print "HE IS STOCHASTICALLY ALIVE!"
                soltype=3
                if savefftanalytics:
                    fftBt = np.absolute((np.fft.fft(bsol_stochastic[:len(bsol_constant)])/n)[range(n/2)])
                    power_fftBt_1yr = np.trapz(np.fabs(fftBt[idx1yr-indx_10delta_frq:idx1yr+indx_10delta_frq]))
                    power_fftBt_0_5yr = np.trapz(np.fabs(fftBt[idx0_5yr-indx_10delta_frq:idx0_5yr+indx_10delta_frq]))
                    power_fftBt_0_25yr = np.trapz(np.fabs(fftBt[idx0_25yr-indx_10delta_frq:idx0_25yr+indx_10delta_frq]))
                    power_fftBt_0_125yr = np.trapz(np.fabs(fftBt[idx0_125yr-indx_10delta_frq:idx0_125yr+indx_10delta_frq]))
                    data['fft2to1ratio_stochastic']=power_fftBt_0_5yr/power_fftBt_1yr
                    data['fft4to1ratio_stochastic']=power_fftBt_0_25yr/power_fftBt_1yr
                    data['fft8to1ratio_stochastic']=power_fftBt_0_125yr/power_fftBt_1yr
                    data['max_around_fft0_5'] = np.amax(fftBt[(data['idx0_5yr']-10*data['indx_10delta_frq']):(data['idx0_5yr']+10*data['indx_10delta_frq'])])
                    first_half_max = np.amax(fftBt[1:(data['idx0_5yr']-10*data['indx_10delta_frq'])])
                    second_half_max = np.amax(fftBt[(data['idx0_5yr']+10*data['indx_10delta_frq']):data['idx1yr']])
                    data['max_besides_fft0_5'] = np.amax([first_half_max,second_half_max])
    data['soltype']=soltype
    return data

def save_check_soltype(p,p_i,chi,chi_i,a,Es,fname,savefftanalytics=False,locker=None,verbose=False):
    if verbose:
        print "Starting p={},chi={}".format(p,chi)
    data = check_soltype(p,chi,a,Es=Es)
    if locker is not None:
        locker.acquire()
    hn.save_p_chi_soltype(fname,p_i,chi_i,p,chi,data,savefftanalytics)
    if locker is not None:
        locker.release()

def apply_async_tlm_p_chi_scan_soltype(pmin,pmax,a,fname,
                                  Nchi=101,dp=0.005,Es=Es,Ps=Ps,
                                  numproc=10,
                                  verbose=False,send_email=None):
    import multiprocessing as mp
    import time
    if send_email is not None:
        import getpass
        pswd = getpass.getpass('Password:')
    if verbose:
        print "Scanning p_min={},p_max={},a={}".format(pmin,pmax,a)
    start = time.time()
    p_range = np.arange(pmin,pmax+dp,dp)
    chi_range = np.linspace(0,1,Nchi)
    hn.setup_p_chi_soltype(fname,Ps,Es,p_range,chi_range)
    from utilities.multiprocess import available_cpu_count
    availble_cpus = int(available_cpu_count() - 2)
    numproc=min(numproc,availble_cpus)
    if verbose:
        print "Using",numproc," processors"
    pool = mp.Pool(processes=numproc)
    m = mp.Manager()
    locker = m.Lock()
    for p_i,p in enumerate(p_range):
        for chi_i,chi in enumerate(chi_range):
            pool.apply_async(save_check_soltype,
                             args = (p,p_i,chi,chi_i,a,Es,fname,locker,verbose))
    pool.close()
    pool.join()
    print "Took ",time.time()-start
    if verbose:
        print "Saved to:",fname
    if send_email is not None:
        try:
            import smtplib
            from socket import gaierror
            server = smtplib.SMTP('smtp.gmail.com', 587)
            server.ehlo()
            server.starttls()
            server.login(send_email, pswd)
            msg = "\r\n".join([
                    "From: {}".format(send_email),
                    "To: {}".format(send_email),
                    "Subject: Apply async calculations finished",
                    "",
                    "Apply async calculations finished and saved to:", fname
                    ])
            server.sendmail(send_email, send_email, msg)
            server.quit()
        except gaierror:
            pass

def save_check_soltype_per_seed(p,p_i,chi,chi_i,seed,seed_i,
                                a,Es,fname,savefftanalytics=False,
                                locker=None,verbose=False):
    Es['stochastic'][-1]=seed
    if verbose:
        print "Starting p={},chi={},seed={}".format(p,chi,seed),' Stochastic:',Es['stochastic']
    data = check_soltype(p,chi,a,savefftanalytics=savefftanalytics,Es=Es)
    if verbose:
        print "For p={},chi={},seed={}".format(p,chi,seed)," soltype=",data['soltype']
    if locker is not None:
        locker.acquire()
    hn.save_p_chi_soltype_per_seed(fname,p_i,chi_i,seed_i,data,savefftanalytics)
    if locker is not None:
        locker.release()
        
def apply_async_tlm_p_chi_scan_soltype_many_seeds(pmin,pmax,a,fname,
                                                  Nchi=101,dp=0.005,Es=Es,Ps=Ps,
                                                  numproc=10,numseeds=100,savefftanalytics=False,
                                                  verbose=False,send_email=None):
    import multiprocessing as mp
    import time
    if send_email is not None:
        import getpass
        pswd = getpass.getpass('Password:')
    if verbose:
        print "Scanning p_min={},p_max={},a={}".format(pmin,pmax,a)
    p_range = np.arange(pmin,pmax,dp)
    chi_range = np.linspace(0,1,Nchi)
    seed_range = np.arange(100,100+numseeds,1)
    hn.setup_p_chi_soltype_many_seeds(fname,Ps,Es,p_range,chi_range,seed_range,savefftanalytics)
    from utilities.multiprocess import available_cpu_count
    availble_cpus = int(available_cpu_count() - 2)
    numproc=min(numproc,availble_cpus)
    if verbose:
        print "Using",numproc," processors"
    pool = mp.Pool(processes=numproc)
    m = mp.Manager()
    locker = m.Lock()
    start = time.time()
    for p_i,p in enumerate(p_range):
        for chi_i,chi in enumerate(chi_range):
            for seed_i,seed in enumerate(seed_range):
#                save_check_soltype_per_seed(p,p_i,chi,chi_i,seed,seed_i,a,Es,
#                                            fname,savefftanalytics,None,verbose)
                pool.apply_async(save_check_soltype_per_seed,
                                 args = (p,p_i,chi,chi_i,seed,seed_i,a,Es,
                                         fname,savefftanalytics,locker,verbose))
    pool.close()
    pool.join()
    finish=time.time()
    if verbose:
        print "Saved to:",fname
        print "Took ",finish-start
    if send_email is not None:
        try:
            import smtplib
            from socket import gaierror
            server = smtplib.SMTP('smtp.gmail.com', 587)
            server.ehlo()
            server.starttls()
            server.login(send_email, pswd)
            msg = "\r\n".join([
                    "From: {}".format(send_email),
                    "To: {}".format(send_email),
                    "Subject: Apply async calculations finished",
                    "",
                    "Apply async calculations finished and saved to:{}. Took {} secs.".format(fname,str(finish-start))
                    ])
            server.sendmail(send_email, send_email, msg)
            server.quit()
        except gaierror:
            pass

def save_analyze_time_series(p,chi,a,i,
                             fname,
                             locker=None,verbose=False):
    if verbose:
        print "Starting p={},chi={}".format(p,chi)
    data = analyze_time_series(p,chi,a,
                               average_stochastic_time_series=100)
    fname=fname+"_a{:2.1f}c{:2.1f}p{:3.2f}".format(a,chi,p).replace(".","_")
    Tcollapse[i]=data['collapse_time']
    if locker is not None:
        locker.acquire()
    dd.save(fname+".hdf5",data,compression='blosc')
    if locker is not None:
        locker.release()


def init(Tcollapse_array):
    global Tcollapse
    Tcollapse=Tcollapse_array

def apply_async_tlm_calc_period_of_collapse(pmin,pmax,chi,a,
                                            dp,fname,numproc=10,
                                            verbose=False,send_email=None):
    print "Calc period of collapse"
    import multiprocessing as mp
    from multiprocessing import Array
    import time
    if send_email is not None:
        import getpass
        pswd = getpass.getpass('Password:')
    if verbose:
        print "Scanning p_min={},p_max={},a={}".format(pmin,pmax,a)
    start = time.time()
    p_range = np.arange(pmin,pmax+dp,dp)
    Tcollapse = Array('d', np.zeros_like(p_range))
    from utilities.multiprocess import available_cpu_count
    availble_cpus = int(available_cpu_count() - 2)
    numproc=min(numproc,availble_cpus)
    if verbose:
        print "Using",numproc," processors"
    pool = mp.Pool(processes=numproc,initializer=init,initargs=(Tcollapse,))
    m = mp.Manager()
    locker = m.Lock()
    for i,p in enumerate(p_range):
        pool.apply_async(save_analyze_time_series,
                         args=(p,chi,a,i,fname,locker,verbose))
    pool.close()
    pool.join()
    print "Took ",time.time()-start
    data = np.array([p_range,np.frombuffer(Tcollapse.get_obj())]).T
    np.savetxt(fname+"_Tcollapse_to_p.dat",data)
    if verbose:
        print "Saved to:",fname
    if send_email is not None:
        try:
            import smtplib
            from socket import gaierror
            server = smtplib.SMTP('smtp.gmail.com', 587)
            server.ehlo()
            server.starttls()
            server.login(send_email, pswd)
            msg = "\r\n".join([
                    "From: {}".format(send_email),
                    "To: {}".format(send_email),
                    "Subject: Apply async calculations finished",
                    "",
                    "Apply async calculations finished and saved to:", fname
                    ])
            server.sendmail(send_email, send_email, msg)
            server.quit()
        except gaierror:
            pass

def main(args):
    if args.analyze_time_series:
        analyze_time_series(args.p,args.chi,args.dist_seed,args.dist_scale,args.a,
                            dist_type=args.dist_type,
                            Es=Es,Ps=args.Ps,
                            savefile=args.fname,timer=args.verbose)
    elif args.scan_p_chi_grid_soltype:
        Es['stochastic']=(args.dist_type,args.scale,args.dist_seed)
        apply_async_tlm_p_chi_scan_soltype(args.pmin,args.pmax,args.a,args.fname,
                                           args.Nchi,args.dp,
                                           Es,args.Ps,args.numproc,
                                           args.verbose,args.send_email)
    elif args.scan_p_chi_grid_soltype_many_seeds:
        Es['stochastic']=[args.dist_type,args.scale,args.dist_seed]
        print Es['stochastic']
        apply_async_tlm_p_chi_scan_soltype_many_seeds(args.pmin,args.pmax,args.a,args.fname,
                                                      args.Nchi,args.dp,
                                                      Es,args.Ps,
                                                      args.numproc,args.average_stochastic_time_series,
                                                      args.savefftanalytics,
                                                      args.verbose,args.send_email)
    elif args.scan_p_chi_grid_analyze:
        apply_async_tlm_p_chi_scan_analyze(args.pmin,args.pmax,args.a,
                                           args.dist_type,args.scale,
                                           args.fname,
                                           args.Nchi,args.dp,
                                           Es,args.Ps,
                                           args.average_stochastic_time_series,
                                           args.numproc,args.savehdf,
                                           args.verbose,args.send_email,args.above_diagonal)
    elif args.scan_period_of_collapse:
        Es['stochastic']=(args.dist_type,args.scale,100)
        apply_async_tlm_calc_period_of_collapse(args.pmin,args.pmax,args.chi,args.a,
                                                args.dp,args.fname,args.numproc,
                                                args.verbose,args.send_email)
    return 0

def add_parser_arguments(parser):
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=False,
                        help="Turn on debuging messages")
    parser.add_argument("--scan_p_chi_grid_soltype",
                        action="store_true",
                        dest="scan_p_chi_grid_soltype",
                        default=False,
                        help="Apply async scan of p to chi grid for soltype only")
    parser.add_argument("--scan_p_chi_grid_soltype_many_seeds",
                        action="store_true",
                        dest="scan_p_chi_grid_soltype_many_seeds",
                        default=False,
                        help="Apply async scan of p to chi grid for soltype only for many seeds")
    parser.add_argument("--savefftanalytics",
                        action="store_true",
                        dest="savefftanalytics",
                        default=False,
                        help="Save savefftanalytics in check_soltype")
    parser.add_argument("--analyze_time_series",
                        action="store_true",
                        dest="analyze_time_series",
                        default=False,
                        help="Run analyze_time_series")
    parser.add_argument("--scan_p_chi_grid_analyze",
                        action="store_true",
                        dest="scan_p_chi_grid_analyze",
                        default=False,
                        help="Apply async scan of p to chi grid for full analysis")
    parser.add_argument("--scan_period_of_collapse",
                        action="store_true",
                        dest="scan_period_of_collapse",
                        default=False,
                        help="Apply async scan of scan_peri--send_emailod_of_collapse")
    parser.add_argument("--savehdf",
                        action="store_true",
                        dest="savehdf",
                        default=False,
                        help="Save final results to a dictionary file using deepdish")
    parser.add_argument("-f", "--fname",
                        type=str, nargs='?',
                        dest="fname",
                        default="tlm_p_chi_scan",
                        help="Save tlm_p_chi_scan in fname")
    parser.add_argument("--Ps",
                        type=str, nargs='?',
                        dest="Ps",
                        default="auto/tlm_set25.hdf5",
                        help="Parameters file name")
    parser.add_argument('--above_diagonal',
                        nargs='+',
                        dest="above_diagonal",
                        default=None,
                        help='Apply async above a certain diagonal')
    parser.add_argument("--dist_type",
                        type=str, nargs='?',
                        dest="dist_type",
                        default="norm_around_p",
                        help="Choose the stochastic probability distribution type")
    parser.add_argument("--dist_seed",
                        type=float, nargs='?',
                        dest="dist_seed",
                        default=100,
                        help="Choose the initial seed")
    parser.add_argument("-s", "--scale",
                        type=float,
                        dest="scale",
                        default=0.1,
                        help="Choose the scale of the probability distribution")
    parser.add_argument("--send_email",
                        type=str, nargs='?',
                        dest="send_email",
                        default=None,
                        help="Password to send email notice")
    parser.add_argument('--pmin',
                        dest='pmin',
                        type=float,
                        default=1.5,
                        help='Minimal Precipitation')
    parser.add_argument('--pmax',
                        dest='pmax',
                        type=float,
                        default=3.5,
                        help='Maximal Precipitation')
    parser.add_argument('-a',
                        dest='a',
                        type=float,
                        default=1.00,
                        help='Seasonality strength')
    parser.add_argument('-p',
                        dest='p',
                        type=float,
                        default=1.00,
                        help='Precipitation')
    parser.add_argument('-c','--chi',
                        dest='chi',
                        type=float,
                        default=0.00,
                        help='Trade off parameter')
    parser.add_argument('--numproc',
                        dest='numproc',
                        type=int,
                        default=10,
                        help='Max number of processors')
    parser.add_argument('--Nchi',
                        dest='Nchi',
                        type=int,
                        default=6,
                        help='Number of points in the chi trade off axis')
    parser.add_argument('--average_stochastic_time_series',
                        dest='average_stochastic_time_series',
                        type=int,
                        default=0,
                        help='Number of runs for averaging fftBt')
    parser.add_argument('--dp',
                        dest='dp',
                        type=float,
                        default=0.01,
                        help='Delta in in the precipitation axis')

    return parser


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    main(parser.parse_args())
