from numpy import pi
dimpar = {'K'   : 1.0,
          'K2'  : 0.4,
          'e'   : 1.0,
          'r'   : 2.0*pi,
          'd'   : 2.0*pi,
          'c'   : 4.0*pi,
          'b'   : 0.3,
          'a'   : 0.0,
          'omegaf'   : 1.0,
         }

def update_par():
    par=dimpar.copy()
    return par


def savecsvdimpar(fname,dimpar):
    import csv
    with open(fname+'.csv', 'wb') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['Parameter','value'])
        writer.writeheader()
        for k in dimpar:
            writer.writerow({'Parameter': k , 'value': dimpar[k]})
    
def saveParmSet(fname,par,text=None,saveroot=False):
    from deepdish.io import save
    if saveroot:
        from tlmModel import tlmModel,Es_normal
        Es_normal['verbose']=False
        m = tlmModel(Es=Es_normal,Ps=par,Vs=None)
        t,result = m.ode_integrate([0.9,0.3,0.3])
        b,s1,s2 = result[-1]
        par['b']=b
        par['s1']=s1
        par['s2']=s2
        print b,s1,s2
    if text is not None:
        par['text']=text
    save("./auto/"+fname+".hdf5",par)
def loadParmSet(fname):
    from deepdish.io import load
    return load(fname)

if __name__ == '__main__':
    par=update_par()
    print "Nondimensional:",
    print par
#    print "conv P=",par['conv_P']
    saveParmSet('rm_set1',par,saveroot=False)
#    import numpy as np
#    p=par
#    a = np.array([p['lamb_max'],p['lamb_min'],p['eta'],
#                  p['p'],p['nu'],p['rho'],p['kappa'],
#                  p['c'],p['zeta'],p['gamma'],
#                  p['s_wp'],p['s_fos'],p['s_fc'],
#                  p['s_h'],p['mu_s_max'],
#                  p['chi'],p['beta'],p['a'],
#                  p['omegaf'],p['delta_s']])
#    np.savetxt("./auto/tlm_parameters.txt",a.T)
##    savemat("./p2p/b2s2_parameters.mat",p)
    
# tlm_set23 is modification of set18 to try to push the bifurcations to lower 
# mm values