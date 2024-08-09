import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.interpolate as inter
import sys
from matplotlib.pyplot import figure

#mDP=float(sys.argv[1])
#print(mDP)

epsilon=np.logspace(-9,-3,50)
totaldecays=np.zeros(len(epsilon))

def chi_calc(m_mu, E_mu, m_z, Z):
    A=184
    Z=74
    me=0.511/1000 #GeV
    a=111*Z**(-1/3)/me
    d=0.164*A**(-2/3)
    t_a=(1/a)**2
    t_d=d
    t_min = m_z**4/(4*E_mu**2)
    t_max = m_z**2 + m_mu**2
    chi_max = t_d**2/(t_a - t_d)**3*((t_a - t_d)*(t_a + t_min)/(t_max + t_a) + (t_a - t_d)*(t_d + t_min)/(t_max + t_d) + (t_a + t_d + 2*t_min)*np.log((t_max+t_d)/(t_max+t_a)))
    chi_min = t_d**2/(t_a - t_d)**3*((t_a - t_d) + (t_a - t_d) + (t_a + t_d + 2*t_min)*np.log((t_min+t_d)/(t_min+t_a)))
    chi_IWW = Z**2*(chi_max - chi_min)
    return(chi_IWW)

def diff_cross_section(eps, alpha, m_z, E_mu, m_mu, Z, x):
    chi_IWW = chi_calc(m_mu, E_mu, m_z, Z)
    u_max = -m_z**2*(1-x)/x - m_mu**2*x
    u_min = -x*(E_mu*0.1)**2 + u_max #theta_max = 0.1 (E_mu(theta) = m_mu/theta)
    dsig_dx_min = 2*eps**2*alpha**3*chi_IWW*np.sqrt((x**2 - (m_z/E_mu)**2))*(m_mu**2*x*(-2+2*x+x**2)-2*(3-3*x+x**2)*u_min)/(3*x*u_min**2)
    dsig_dx_max = 2*eps**2*alpha**3*chi_IWW*np.sqrt((x**2 - (m_z/E_mu)**2))*(m_mu**2*x*(-2+2*x+x**2)-2*(3-3*x+x**2)*u_max)/(3*x*u_max**2)
    dsig_dx = dsig_dx_max - dsig_dx_min

    return(dsig_dx)


N=1000
Z=74
m_mu = 0.1057 #GeV/c^2
alpha = 1/137
L_target=0.5 #cm
me=0.511/1000 #GeV
mp=938.3/1000 #GeV
E_e = me+0.1 #GeV
E_p = mp+0.016 #GeV
Na=6.02*10**23
A=184

mDP=np.logspace(-2.987,-1,30)
dsig_dx=np.zeros((len(mDP),N))
ratio=np.zeros(N)

for z in range(len(mDP)):
    print(mDP[z])
    f=open('/home/LuisFranca/forDavid/times/BaF2/ndecays_%s_BaF2_eedecay_yemilab.txt' % "{:.2e}".format(mDP[z]),'w')
    for k in range(len(epsilon)):
        E_z = np.logspace(np.log10(mDP[z]), np.log10(E_e), N) #GeV
        totsig=0.0
        for r in range(N):
            ratio[r] = E_z[r]/E_e
            if (r!=0):
                dsig_dx[z][r] = diff_cross_section(epsilon[k], alpha, mDP[z], E_e, me, Z, (ratio[r]+ratio[r-1])/2)
                integr=(dsig_dx[z][r]+dsig_dx[z][r-1])*(ratio[r]-ratio[r-1])/2
                sig=integr*(1.973*10**(-16))**2*10**4*10**24*10**12 #pb
                #totsig=totsig+sig
                xsec=sig*10**(-12)*10**(-24) #cm^2

                Nbremss=6.2E+15*11.3*Na/A*L_target*xsec*60*60*24 #day^-1

                L=1.973E-16*2*(E_z[r]+E_z[r-1])/2/(epsilon[k]**2*alpha*mDP[z]**2*(1+2*me**2/mDP[z]**2)*(1-4*me**2/mDP[z]**2)**0.5) #m
                d=0.5 #m (distance between target and BaF2)
                P=np.exp(-(d/L))-np.exp(-((d+20)/L))
                totaldecays[k]=totaldecays[k]+P*Nbremss

        #print(mDP[z],epsilon[k],totaldecays[k],totsig)
        print("{:.3e}".format(epsilon[k]),file=f,end=" ")
        print("{:.3e}".format(totaldecays[k]),file=f,end="\n")
    f.close()