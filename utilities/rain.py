# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 22:50:02 2016
http://web.mit.edu/1.017/www/lecnotes_03/extras/Poisson/Poisson00.html
@author: ohm
"""

def main(period,rspacing,rmean):
    """
    n=100       % number of days simulated
    rspacing=2  % average days between rain
    rmean=10.   % average rainfall on rainy days
    time=[0:n];
#% call rain generator
#precip=raingen(n,rspacing,rmean);
#% augment precip series for plotting
#precip=[precip,precip(n)];
#% plot precip series
#close all
#figure
#stairs(time,precip)
#title('Simulated Rainfall Series')
#axis([0 time(n+1) 0 1.1*max(precip)])
#xlabel('Time')
#ylabel('Precip., mm/day')
#return
#%
    """
    return 0
    
def poisson_series(ndays, spacing, meanrain):
    """
    #% Raingen generates a synthetic (random) sequence
    #% of daily rainfall events.
    #% Event spacing is Poisson distributed
    #% Event magnitude is exponentially distributed
    #% ndays is total number of days simulated
    #% spacing is mean spacing between rain days
    #% meanrain is mean daily rainfall
    #% rain is array of generated daily rainfall  values
    """

#% raindays is array of rain day numbers
#te=1;
#% choose first rainday (must be >0) from Poisson generator
#intvl=poissrnd(spacing,1,1);
#rday=1+intvl;
#% generate other raindays by adding number of days until
#% next rain (obtained from Poisson RNG)
#% continue until final day has been passed
#while rday<=ndays
#   raindays(te)=rday;
#   intvl=poissrnd(spacing,1,1)
#   rday=rday+intvl;
#   te=te+1;
#end
#% initialize output array (default rainfall is zero)
#nrain=te-1;
#rain=zeros(1,ndays);
#% generate exponentially distributed rainfall
#% depths at raindays
#rain(raindays)=exprnd(meanrain,1,nrain);
#return
if __name__ == '__main__':
    main()