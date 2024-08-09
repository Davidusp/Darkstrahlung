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

# Calculate cross section considering Earth layers
@numba.jit(nopython=True, nogil=True)
def chi_calc(m_mu, E_mu, m_z, Z, A):
    me=0.511/1000 #GeV
    a=111*Z**(-1/3)/me
    d=0.164*A**(-2/3)
    t_a=(1/a)**2
    t_d=d
    #t_a = (2e-5)**2 #GeV^2
    #t_d = (67e-3)**2 #GeV^2
    t_min = m_z**4/(4*E_mu**2)
    t_max = m_z**2 + m_mu**2
    chi_max = t_d**2/(t_a - t_d)**3*((t_a - t_d)*(t_a + t_min)/(t_max + t_a) + (t_a - t_d)*(t_d + t_min)/(t_max + t_d) + (t_a + t_d + 2*t_min)*np.log((t_max+t_d)/(t_max+t_a)))
    chi_min = t_d**2/(t_a - t_d)**3*((t_a - t_d) + (t_a - t_d) + (t_a + t_d + 2*t_min)*np.log((t_min+t_d)/(t_min+t_a)))
    chi_IWW = Z**2*(chi_max - chi_min)
    return(chi_IWW)

@numba.jit(nopython=True, nogil=True)
def diff_cross_section(eps, alpha, m_z, E_mu, m_mu, Z, A):
    '''
        Calculate the differential cross section for DP production by muon Bremsstrahlung
        using the IWW approximation
    '''
    m_z = m_z/1000
    dsig_dx = []

    for i in range(0, len(E_mu)):
        E_z = np.linspace(m_z, E_mu[i], 1000)
        chi_IWW = chi_calc(m_mu, E_mu[i], m_z, Z, A)
        x = E_z/E_mu[i]
        u_max = -m_z**2*(1-x)/x - m_mu**2*x
        u_min = -x*(E_mu[i]*0.1)**2 + u_max #theta_max = 0.1 (E_mu(theta) = m_mu/theta)
        #print(u_max, u_min)
        dsig_dx_min = 2*eps**2*alpha**3*chi_IWW*np.sqrt((x**2 - (m_z/E_mu[i])**2))*(m_mu**2*x*(-2+2*x+x**2)-2*(3-3*x+x**2)*u_min)/(3*x*u_min**2)
        dsig_dx_max = 2*eps**2*alpha**3*chi_IWW*np.sqrt((x**2 - (m_z/E_mu[i])**2))*(m_mu**2*x*(-2+2*x+x**2)-2*(3-3*x+x**2)*u_max)/(3*x*u_max**2)
        dsig_dx.append(list(dsig_dx_max - dsig_dx_min))

    return(np.array(dsig_dx))

def main():

    m_mu = 0.105658     #GeV  #muon mass
    Z = 11              #Rock atomic number
    A = 22              #Rock mass number
    m_z = 1e01          #DP mass (MeV)
    epsilon = 1e-05     #coupling
    alpha = 1/137

    E_mu = np.array([1e-01, 1, 10, 100, 1000]) #muon energy (GeV)
    
    dsig_dx = diff_cross_section(epsilon, alpha, m_z, E_mu, m_mu, Z, A)

    x = np.linspace(m_z/(E_mu[3]*1e03), 1, 1000)
 
    E1 = np.linspace(m_z/(E_mu[0]*1e03), 1, 1000)
    E2 = np.linspace(m_z/(E_mu[1]*1e03), 1, 1000)
    E3 = np.linspace(m_z/(E_mu[2]*1e03), 1, 1000)
    E4 = np.linspace(m_z/(E_mu[3]*1e03), 1, 1000)
    E5 = np.linspace(m_z/(E_mu[4]*1e03), 1, 1000)

    ds_dx1 = inter.interp1d(E1, dsig_dx[0], kind="cubic", fill_value="extrapolate")
    ds_dx2 = inter.interp1d(E2, dsig_dx[1], kind="cubic", fill_value="extrapolate")
    ds_dx3 = inter.interp1d(E3, dsig_dx[2], kind="cubic", fill_value="extrapolate")
    ds_dx4 = inter.interp1d(E4, dsig_dx[3], kind="cubic", fill_value="extrapolate")
    ds_dx5 = inter.interp1d(E5, dsig_dx[4], kind="cubic", fill_value="extrapolate")

    #Plotting diferential cross-sections
    f1=plt.figure()
    plt.title("Differential bremsstrahlung cross-section - $m_{Z'}$ = 10 MeV")
    plt.plot(x, ds_dx1(x)/epsilon**2, color='blue', linestyle='-', linewidth=1.5, label='$E_{\mu} = 100 MeV$')
    plt.plot(x, ds_dx2(x)/epsilon**2, color='black', linestyle='-.', linewidth=1.5, label='$E_{\mu} = 1 GeV$')
    plt.plot(x, ds_dx3(x)/epsilon**2, color='green', linestyle='--', linewidth=1.5, label='$E_{\mu} = 10 GeV$')
    plt.plot(x, ds_dx4(x)/epsilon**2, color='red', linestyle='-.', linewidth=1.5, label='$E_{\mu} = 100 GeV$')
    plt.plot(x, ds_dx5(x)/epsilon**2, color='orange', linestyle=':', linewidth=1.5, label='$E_{\mu} = 1 TeV$')
    plt.tick_params(axis='both', labelsize=12)  # Tick labels fontsize
    plt.xlabel("$x$", fontsize=13)
    plt.ylabel("$d\sigma/(\epsilon^2 dx) [GeV^{-2}]$", fontsize=13)
    plt.yscale('log')
    plt.xscale('linear')
    #plt.ylim(1E+2, 5e+9)
    #plt.xlim(1E-2,1)
    #plt.legend(loc = 'upper right', fontsize="12")
    plt.grid(True)
    plt.savefig("/home/davidc/Documents/Master's Analisys/Parameter space/Codes/dsig_dx_m=%.8f_k=%.11f.png" %(m_z, epsilon))
    #plt.show()
    plt.close()

main()