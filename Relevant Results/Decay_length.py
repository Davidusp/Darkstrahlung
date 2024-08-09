import math
from math import isnan, nan
from os import O_SYNC, kill
from tracemalloc import stop
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.interpolate as inter
import sys
import numba

def CalcpairsL(mDP,kappa, E_z, alpha, m_e):

    gamma=kappa**2*alpha*mDP/3*(1+2*m_e**2/mDP**2)*(1-4*m_e**2/mDP**2)**0.5 #decay rate
    Ldec = np.sqrt(E_z**2/mDP**2 - 1)/gamma  #decay length (MeV⁻¹)

    return(Ldec)

def Calc3gammaL(mDP,kappa, E_z, alpha, m_e):

    gamma = kappa**2*alpha**4/(2**7*3**6*5**2*np.pi**3)*mDP**9/m_e**8*(17/5+67/42*mDP**2/m_e**2+128941/246960*mDP**4/m_e**4)
    Ldec = np.sqrt(E_z**2/mDP**2 - 1)/gamma  #decay length (MeV⁻¹)

    return(Ldec)

def main():

    m_e = 0.511 #electron mass (MeV)
    alpha = 1/137

    mDP1 = np.logspace(np.log10(2*m_e+1e-05*m_e), np.log10(10), 100)
    kappa1 = 1e-05
    mDP2 = np.logspace(np.log10(1e-2), np.log10(1), 100)
    kappa2 = 1e-01

    E_z1 = 100
    E_z2 = 2

    Lpairs = CalcpairsL(mDP1, kappa1, E_z1, alpha, m_e)*1.98e-13 #(m)
    L3gamma = Calc3gammaL(mDP2, kappa2, E_z2, alpha, m_e)*1.98e-13 #(m)

    f1=plt.figure()
    plt.title("Decay Length ($e^{+}e^{-}$), $\epsilon=10^{-5}$")
    plt.plot(mDP1,Lpairs,c='blue')
    plt.xlabel("$m_{Z'} [MeV]$", fontsize=11)
    plt.ylabel("Decay Length [m]", fontsize=12)
    plt.yscale('log')
    plt.xscale('log')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=11)
    #plt.legend(loc = 'lower right', fontsize="9")
    plt.grid(True)
    plt.savefig("/home/davidc/Documents/Master's Analisys/Parameter space/Codes/DecayLpairs_E=%f_k=%f.png" %(E_z1, kappa1))
    #plt.show()
    plt.close()

    f1=plt.figure()
    plt.title("Decay Length ($3\gamma$), $\epsilon=10^{-1}$")
    plt.plot(mDP2,L3gamma,c='blue')
    plt.xlabel("$m_{Z'} [MeV]$", fontsize=11)
    plt.ylabel("Decay Length [m]", fontsize=12)
    plt.yscale('log')
    plt.xscale('log')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=11)
    #plt.legend(loc = 'lower right', fontsize="9")
    plt.grid(True)
    plt.savefig("/home/davidc/Documents/Master's Analisys/Parameter space/Codes/DecayL3g_E=%f_k=%f.png" %(E_z2, kappa2))
    #plt.show()
    plt.close()

main()