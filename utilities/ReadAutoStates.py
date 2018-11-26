# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 11:01:20 2016

@author: ohm
"""
import numpy as np
from scipy.interpolate import interp1d
def ReadAutoStates(fname,Ps,Es={},pnts=None,varargin=None):
    if 'badtext' in Es:
        Es['badtext']=0
    #Initialization
    if pnts is None:
        pnts = 0
    abit = 0.001       # a very small value, used for interpolation 
    linelen = 16       # the number of parms in an info line (before each state)
    vrnum = Ps['Vnum'] # how many variables to read per row
    fsize = Ps['Nx']   # how many points to interpolate into
    maxpnt = np.max(pnts) # how many different states to try and read
    
    # read file and read the info-line of the first state
    with open(fname, 'r') as fin_org:
        fin = fin_org.read().split()
    infoline = [int(fin.pop(0)) for i in range(linelen)]#map(int,fin[:linelen])
    rsz = infoline[6] # number of points (rows) saved per state
    wid = infoline[7] # number of variables (columns)
    totrowlen = infoline[8]+1 # how many rows per state in total - unused
    extra = infoline[11] # number of extra numbers at the end 
    lines = infoline[14] # number of lines this takes
    lenfin = len(fin)
    num_throw = rsz*(wid-1)+extra
    print num_throw
    print "first info line", infoline
#        tmpst = np.zeros(wid,rsz) #  this will hold the data for each state we read
    
    states = []
    indx   = []
    ind = 0 
    while (len(infoline)==linelen) and ((pnts==0) or (ind<maxpnt)):
        indx.append(infoline[1]) # Add this state's index
        # read data of state
        state = []
        for i in range(rsz):
            state.append([float(fin.pop(0)) for j in range(wid)])
        state = np.array(state)
        Xs = np.hstack((0-abit,state[:,0],1+abit))
        Ys = np.vstack((state[0,1:vrnum+1],state[:,1:vrnum+1],state[-1,1:vrnum+1]))
        states.append(np.array(interp1d(Xs,Ys, kind='cubic',axis=0)(np.linspace(0,1,fsize))).T)
        # read through data that comes at the end of each state, that we do not use
        for i in range(num_throw):
            temp = fin.pop(0)
#            print i, temp
        infoline = [int(fin.pop(0)) for i in range(linelen)]
        print infoline
        ind+=1
            
#        if pnts!=0:
#            states = states[pnts]
    return states


#function [states, indx]=ReadAutoStates(filename,Ps,Es,pnts,varargin)
#% Read AUTO states file (usually s.name or fort.8)
#% points=ReadAutoStates(filename,Ps,Es,pnts)
#% Returns a presentation of a set of states, identified by pnts
#
#% Update online if necessary
#[~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});
#
#if(~isfield(Es,'badtext')) % fix AUTO's bad-text problem (for very small/big values)
#   Es.badtext=0;
#end;
#
#% Initialization
#if(nargin<4) 
#	pnts=0;
#end;
#abit    = 0.001;        % a very small value, used for interpolation
#linelen = 16;           % the number of parms in an info line (before each state)
#vrnum   = Ps.Vnum;      % how many variables to read per row
#fsize   = Ps.Nx;        % how many points to interpolate into
#maxpnt  = max(pnts);    % how many different states to try and read
#
#% read file and read the info-line of the first state
#fin   = fopen(filename,'r');
#infoline  = fscanf(fin,'%f',linelen);
#
#% get some general parameter from this infoline
#rsz   = infoline(7);    % number of points (rows) saved per state
#wid   = infoline(8);    % number of variables (columns)
#totrowlen  = infoline(9)+1;  % how many rows per state in total - unused
#extra = infoline(5)*2+infoline(12);    % number of extra numbers at the end 
#lines = infoline(14)+infoline(15);     % number of lines this takes
#
#tmpst = zeros(wid,rsz); % this will hold the data for each state we read
#
#ind    = 0;
#states = [];
#indx   = [];
#
#
#
#while((length(infoline)==linelen) && ((pnts(1)==0)||(ind<maxpnt)))
#	indx = [indx infoline(2)]; % Add this state's index
#    % read data of state
#    if(Es.badtext==0)
#    	[tmpst(:),~] = fscanf(fin,'%f',rsz*wid);    % Faster better method
#    else
#        tmpst = ReadSingleState(fin,rsz,wid)';  % Much slower, but works with bad files
#    end;
#	% read things for interpolation, and then do it
#	Xs = [0-abit  tmpst(1,:)  1+abit ]';
#	Ys = [tmpst(2:vrnum+1,1) tmpst(2:vrnum+1,:) tmpst(2:vrnum+1,end)]';
#	states = cat(3,states,interp1(Xs,Ys,(0:(fsize-1))'/fsize));
#	
#
#    % read through data that comes at the end of each state, that we do not use
#    if(Es.badtext==0)
#    	 fscanf(fin,'%f',rsz*(wid-1)+extra);
#    else
#        for ii=1:(rsz+lines) fgetl(fin); end;
#    end;
#	
#    % read next infoline of new state
#    infoline  = fscanf(fin,'%f',linelen);
#
#	ind = ind+1;
#end;
#% if not all points were asked for, then return only the specific list
#if(pnts~=0)
#	states = states(:,:,pnts);
#end;
#
#fclose(fin);
#
#end