import numpy as np
import sys
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.interpolate as inter
from scipy.integrate import quad
import numba
from numba import jit

#@numba.jit(nopython=True, nogil=True)
def CalcAbs(mDP, kappas):
    '''
    Function to estimate the number of absorptions occuring
    for each coupling parameter and DP mass
    '''
    q = 0.30282212088 #Elementary charge in natural units
    alpha = q**2/(4*np.pi) #Fine structure constant
    alphaline = (q*kappas)**2/(4*np.pi)
    Na=6.02*10**23    #Avogadro number

    # with open('/home/dfreitas/DPnum_Comp/Absorption/PeXS_E vs Sig.txt', 'r') as archive:
    #         arq = archive.read() 

    with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/PeXSNaI.txt', 'r') as archive:
            arq = archive.read()

    lines = arq.split("\n")
    lines = lines[:-1]
    
    for row in range(0, len(lines)):
        lines[row] = lines[row].split(" ")
        lines[row] = [float(x) for x in lines[row]]

    sigAbs = np.zeros(len(lines))
    E_k = np.zeros(len(lines))

    for l in range(0, len(lines)):
        E_k[l] = lines[l][0]
        if mDP < E_k[l]:
            sigAbs[l] = (lines[l][1]/np.sqrt(1-mDP**2/E_k[l]**2))*(alphaline/alpha)*(149.89/Na) *(1e24) #barns/atom

    #CSabsfunc = inter.interp1d(E_k, sigAbs, kind = 'linear', fill_value="extrapolate")
    CSabsfunc = inter.interp1d(E_k, sigAbs, kind = 'linear',fill_value="extrapolate")

    E_k = np.logspace(np.log10(mDP+0.001*mDP),np.log10(100), len(E_k))
    
    sigAbs = CSabsfunc(E_k)
    
    return(sigAbs, E_k)


def main():

    #N = 17000
    kappa = 0.0001  #knitic mixing parameter
    #v=220*17000/c
    #print(v)
    #m_phi = float(sys.argv[1])  #MeV
    m_phi = 1.000000 #MeV

    (sig, E_k) = CalcAbs(m_phi,kappa)
    #E_k = np.logspace(np.log10(m_phi),2,len(sig))  #MeV

    #sig = Calc_sig(dsig_dcos, sig, ppabs, pp, p, qabs, q, kabs, costheta, E_k, m_e, m_phi, kappa, alpha)

    # with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/CrossSec/CSxE_m=%f_k=%f_NaI.txt' %(m_phi, kappa), 'w') as archive:
    #     for i in range(0, len(E_k)):
    #         archive.write('%.8f   %.16f\n' %(E_k[i], sig[i]))

    plt.title('Absorption Cross-section')
    plt.plot(E_k, sig, c='blue', linewidth=2)
    plt.xlabel("$E_{k}[MeV]$",fontsize=12)
    plt.ylabel("$\sigma~[barns/atom]$",fontsize=12)
    plt.yscale('log')
    plt.xscale('linear')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.ylim(1E+2, 5e+9)
    #plt.xlim(1E-2,1)
    plt.legend(loc = 'upper right')
    plt.grid(True)
    plt.show()
    #plt.savefig("/home/dfreitas/Results/CS_%dE_%dcos_k=%f_mphi=%f.png" %(len(E_k), len(q), kappa, m_phi))

main()
