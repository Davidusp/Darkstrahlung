import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.interpolate as inter
from scipy.integrate import quad
import numba
from numba import jit

def calc_matrix_el(E_pp, qabs, pp, p, q, m_e, m_phi):

    '''
       This function computes the interaction amplitude matrix
       element for a 2 to 2 proccess with a vectorial boson (A_v).
    '''

    #Mandelstam variables
    s = -2*(-E_pp*qabs+np.dot(pp,q))   #-(pp + k)**2 - m_e**2
    u = 2*(-m_e*qabs+np.dot(p,q))     #-(p - k)**2 - m_e**2
    t_2 = 2*(-E_pp*m_e+np.dot(pp,p)) + 2*m_e**2  #-(pp - p)**2
        
    A_v = 4 - 2*(s+u)**2/(s*u) + 4*(m_phi**2 + 2*m_e**2)*(((s+u)/(s*u))**2*m_e**2 - t_2/(s*u))

    return(A_v)

@numba.jit(nopython=True, nogil=True)
def Calc_kq(v, kpre, k, kabs, costheta, qpre, q, qabs, m_e, m_phi):

    for m in range(len(m_phi)):
        print(m)
        i=0
        while (i<1000):
            kpre[0]=np.random.normal(0,m_phi[m]*v)
            kpre[1]=np.random.normal(0,m_phi[m]*v)
            kpre[2]=np.random.normal(0,m_phi[m]*v)
            kpreabs=(kpre[0]**2+kpre[1]**2+kpre[2]**2)**0.5
            if (kpreabs>(m_phi[m]*v-m_phi[m]*v*0.01) and kpreabs<(m_phi[m]*v+m_phi[m]*v*0.01)):
                kabs[m][i]=kpreabs
                k[m][i][0]=kpre[0]
                k[m][i][1]=kpre[1]
                k[m][i][2]=kpre[2]
                i=i+1
    
        costheta[m]=np.linspace(-1,1,1000)

        i=0
        while (i<1000):
            qpre[0]=np.random.normal(0,(m_phi[m]**2+2*m_phi[m]*m_e)/(2*(m_e+m_phi[m]-kabs[m][i]*costheta[m][i])))
            qpre[1]=np.random.normal(0,(m_phi[m]**2+2*m_phi[m]*m_e)/(2*(m_e+m_phi[m]-kabs[m][i]*costheta[m][i])))
            qpre[2]=np.random.normal(0,(m_phi[m]**2+2*m_phi[m]*m_e)/(2*(m_e+m_phi[m]-kabs[m][i]*costheta[m][i])))
            qpreabs=(qpre[0]**2+qpre[1]**2+qpre[2]**2)**0.5
            if (qpreabs>((m_phi[m]**2+2*m_phi[m]*m_e)/(2*(m_e+m_phi[m]-kabs[m][i]*costheta[m][i])) - (m_phi[m]**2+2*m_phi[m]*m_e)/(2*(m_e+m_phi[m]-kabs[m][i]*costheta[m][i]))*0.01) and qpreabs<((m_phi[m]**2+2*m_phi[m]*m_e)/(2*(m_e+m_phi[m]-kabs[m][i]*costheta[m][i])) + (m_phi[m]**2+2*m_phi[m]*m_e)/(2*(m_e+m_phi[m]-kabs[m][i]*costheta[m][i]))*0.01)):
                if (np.dot(k[m][i],qpre)/(qpreabs*kabs[m][i])>(costheta[m][i]-0.02) and np.dot(k[m][i],qpre)/(qpreabs*kabs[m][i])<(costheta[m][i]+0.02)):
                    #print(m,i)
                    qabs[m][i]=qpreabs
                    q[m][i][0]=qpre[0]
                    q[m][i][1]=qpre[1]
                    q[m][i][2]=qpre[2]
                    i=i+1

    return(q, qabs)

@numba.jit(nopython=True, nogil=True)
def Calc_p(m_phi, k, q, pp, ppabs):

    for m in range(len(m_phi)):
        for i in range(1000):
            pp[m][i][0]=k[m][i][0]-q[m][i][0]
            pp[m][i][1]=k[m][i][1]-q[m][i][1]
            pp[m][i][2]=k[m][i][2]-q[m][i][2]
            ppabs[m][i]=(pp[m][i][0]**2+pp[m][i][1]**2+pp[m][i][2]**2)**0.5
            #print(m,np.dot(k[m][i],q[m][i])/(qabs[m][i]*kabs[m][i]),costheta[m][i])

    return(pp, ppabs)

def CalcAbs(kappas):
    '''
    Function to estimate the number of Compton-like scatterings occuring for each angle
    and step inside the rock, and for each coupling parameter and DP mass 
    '''
    q = 0.30282212088 #Elementary charge in natural units
    alpha = 1/137 #q**2/(4*np.pi) #Fine structure constant
    alphaline = (q*kappas)**2/(4*np.pi)
    v = 10**(-3) #speed of DP in CDM model
    Na = 6.022e23 #Avogadro number
    AGe = 72.64 #mass number of Ge
    with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/PeGeXS.txt', 'r') as archive:
            arq = archive.read() 

    # with open('/home/dfreitas/DPnum_Comp/Comp-xs/CSxE_m=%f_k=%s00.txt' %(mDP, epsilon), 'r') as archive:
    #         arq = archive.read()

    lines = arq.split("\n")
    lines = lines[:-1]

    for row in range(0, len(lines)):
        lines[row] = lines[row].split(" ")
        lines[row] = [float(x) for x in lines[row]]

    sigAbs = np.zeros(len(lines))
    E_phot = np.zeros(len(lines))

    for l in range(0, len(lines)):
        E_phot[l] = lines[l][0]
        sigAbs[l] = (lines[l][1]/v)*(alphaline/alpha)*1e24*(AGe/Na)#*(1e24/8.36e21) #barns/atom

    CSabsfunc = inter.interp1d(E_phot, sigAbs, kind = 'linear')
    #CSabsfunc = inter.interp1d(E_phot, sigAbs, kind = 'linear', fill_value="extrapolate")

    return(CSabsfunc)

    #return(E_phot,sigAbs)

def main():

    #N = 2000
    c = 299792458 #m/s
    kappa = 1
    v=220*1000/c
    #print(v)
    alpha = 1/137
    m_e = 0.511  #MeV/c²
    m_phi = np.logspace(np.log10(1e-3),np.log10(1),36) #MeV/c²
    
    #E_k=m_phi

    #k = m_phi*v
    kpre=np.zeros(3)
    k=np.zeros((len(m_phi),1000,3))
    kabs=np.zeros((len(m_phi),1000))
    
    qpre=np.zeros(3)
    q=np.zeros((len(m_phi),1000,3))
    qabs=np.zeros((len(m_phi),1000))

    costheta=np.zeros((len(m_phi),1000))
      
    pp=np.zeros((len(m_phi),1000,3))
    ppabs=np.zeros((len(m_phi),1000))
        
    p=[0.0,0.0,0.0]
    dsig_dcos=np.zeros((len(m_phi),1000))
    sigC=np.zeros(len(m_phi))

    (q, qabs) = Calc_kq(v, kpre, k, kabs, costheta, qpre, q, qabs, m_e, m_phi)
    (pp, ppabs) = Calc_p(m_phi, k, q, pp, ppabs)

    for m in range(len(m_phi)):
        #print(j)
        for i in range(1000): 
            E_pp=(ppabs[m][i]**2+m_e**2)**0.5

            A_v = calc_matrix_el(E_pp, qabs[m][i], pp[m][i], p, q[m][i], m_e, m_phi[m])
        
            dsig_dcos[m][i] = qabs[m][i]*A_v/(m_phi[m] + m_e - kabs[m][i]*costheta[m][i])
        
            #if (i%100==0):
            #    print(m,i,A_v)
    
        func=inter.interp1d(costheta[m],dsig_dcos[m])

        def integrand(x):
            return func(x)
        I = quad(integrand, costheta[m][0], costheta[m][len(costheta[m])-1]) 

        Ne=32
        sigC[m] = Ne*kappa**2*np.pi*alpha**2/(3*m_e*m_phi[m]*v)*I[0]
        sigC[m]=sigC[m]*(1.973*10**-11)**2 #MeV^-2 to cm^2
        sigC[m]=sigC[m]*10**24 #cm^2 to barn

    #(mDPs,sigA) = CalcAbs(kappa)
    CSabsfunc = CalcAbs(kappa)

    plt.title('Cross-section')
    plt.plot(m_phi, sigC, c='blue', linewidth=2, linestyle='--', label='Compton-like process')
    plt.plot(m_phi, CSabsfunc(m_phi), c='blue', linewidth=2, linestyle=':', label='Absorption process')
    plt.plot(m_phi, sigC+CSabsfunc(m_phi), c='blue', linewidth=2, label='Total Cross Section')
    plt.xlabel("$m_{\phi}[MeV]$",fontsize=13)
    plt.ylabel("$\sigma~[barn/atom]$",fontsize=13)
    plt.yscale('log')
    plt.ylim(1E+2, 5e+9)
    plt.xlim(0,1)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.legend(loc = 'upper right')
    plt.grid(True)
    plt.show()

main()