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

def calc_matrix_el(E_pp, qabs, pp, p, q, m_e, m_phi):

    '''
       This function computes the interaction amplitude matrix
       element for a 2 to 2 proccess with a vectorial boson (A_v).
    '''

    #Mandelstam variables
    s = -2*(-E_pp*qabs+np.dot(pp,q))                   #-(pp + k)**2 - m_e**2
    u = 2*(-m_e*qabs+np.dot(p,q))                      #-(p - k)**2 - m_e**2
    t_2 = 2*(-E_pp*m_e+np.dot(pp,p)) + 2*m_e**2        #-(pp - p)**2
        
    A_v = 4 - 2*(s+u)**2/(s*u) + 4*(m_phi**2 + 2*m_e**2)*(((s+u)/(s*u))**2*m_e**2 - t_2/(s*u))
    return(A_v)

@numba.jit(nopython=True, nogil=True)
def Calc_kq(E_k, kpre, k, kabs, costheta, qpre, q, qabs, m_e, m_phi):
    for m in range(len(E_k)):
        print(m)
        i=0
        eps_DP = (E_k[m]**2-m_phi**2)**0.5*0.0007
        while (i<17000):
            kpre[0]=np.random.normal(0,((E_k[m]**2-m_phi**2)**0.5))
            kpre[1]=np.random.normal(0,((E_k[m]**2-m_phi**2)**0.5))
            kpre[2]=np.random.normal(0,((E_k[m]**2-m_phi**2)**0.5))
            kpreabs=(kpre[0]**2+kpre[1]**2+kpre[2]**2)**0.5
            if (kpreabs<(E_k[m]**2-m_phi**2)**0.5 and kpreabs>((E_k[m]**2-m_phi**2)**0.5-eps_DP)):
                kabs[m][i]=kpreabs
                k[m][i][0]=kpre[0]
                k[m][i][1]=kpre[1]
                k[m][i][2]=kpre[2]
                i=i+1
    
        costheta[m]=np.linspace(-1,1,17000)

        i=0
        while (i<17000):
            qpre[0]=np.random.normal(0,(m_phi**2+2*E_k[m]*m_e)/(2*(m_e+E_k[m]-kabs[m][i]*costheta[m][i])))
            qpre[1]=np.random.normal(0,(m_phi**2+2*E_k[m]*m_e)/(2*(m_e+E_k[m]-kabs[m][i]*costheta[m][i])))
            qpre[2]=np.random.normal(0,(m_phi**2+2*E_k[m]*m_e)/(2*(m_e+E_k[m]-kabs[m][i]*costheta[m][i])))
            qpreabs=(qpre[0]**2+qpre[1]**2+qpre[2]**2)**0.5

            eps_q = (m_phi**2+2*E_k[m]*m_e)/(2*(m_e+E_k[m]-kabs[m][i]*costheta[m][i]))*0.007
            eps_cos = 0.007
            if (qpreabs<(m_phi**2+2*E_k[m]*m_e)/(2*(m_e+E_k[m]-kabs[m][i]*costheta[m][i])) and qpreabs>((m_phi**2+2*E_k[m]*m_e)/(2*(m_e+E_k[m]-kabs[m][i]*costheta[m][i])) - eps_q)):
                if (np.dot(k[m][i],qpre)/(qpreabs*kabs[m][i])>(costheta[m][i]-eps_cos) and np.dot(k[m][i],qpre)/(qpreabs*kabs[m][i])<(costheta[m][i]+eps_cos)):
                    #print(m,i)
                    qabs[m][i]=qpreabs
                    q[m][i][0]=qpre[0]
                    q[m][i][1]=qpre[1]
                    q[m][i][2]=qpre[2]
                    i=i+1
    return(q, qabs)

@numba.jit(nopython=True, nogil=True)
def Calc_p(E_k, k, q, pp, ppabs):

    for m in range(len(E_k)):
        for i in range(17000):
            pp[m][i][0]=k[m][i][0]-q[m][i][0]
            pp[m][i][1]=k[m][i][1]-q[m][i][1]
            pp[m][i][2]=k[m][i][2]-q[m][i][2]
            ppabs[m][i]=(pp[m][i][0]**2+pp[m][i][1]**2+pp[m][i][2]**2)**0.5
            #print(m,np.dot(k[m][i],q[m][i])/(qabs[m][i]*kabs[m][i]),costheta[m][i])
    return(pp, ppabs)

def main():

    #N = 17000
    #c = 299792458 #m/s
    kappa = 1e-4  #knitic mixing parameter
    #v=220*17000/c
    #print(v)
    alpha = 1/137 #fine structure constant
    m_e = 0.511  #MeV/c²
    #m_phi = float(sys.argv[1])  #MeV
    m_phi = 0.010985
    E_k = np.logspace(np.log10(m_phi+0.001*m_phi),np.log10(100),80) # DP energy (MeV/c²)

    #k = sqrt(E_k²-m_phi²)
    kpre = np.zeros(3)
    k = np.zeros((len(E_k),17000,3))
    kabs = np.zeros((len(E_k),17000))

    qpre = np.zeros(3)
    q = np.zeros((len(E_k),17000,3))
    qabs = np.zeros((len(E_k),17000))

    pp=np.zeros((len(E_k),17000,3))
    ppabs=np.zeros((len(E_k),17000))

    costheta = np.zeros((len(E_k),17000))

    p=[0.0,0.0,0.0] 
    dsig_dcos=np.zeros((len(E_k),17000))
    sig=np.zeros(len(E_k))

    (q, qabs) = Calc_kq(E_k, kpre, k, kabs, costheta, qpre, q, qabs, m_e, m_phi)
    (pp, ppabs) = Calc_p(E_k, k, q, pp, ppabs)

    for m in range(len(E_k)):
        #print(j)
        for i in range(17000): 
            E_pp=(ppabs[m][i]**2+m_e**2)**0.5
        
            A_v = calc_matrix_el(E_pp, qabs[m][i], pp[m][i], p, q[m][i], m_e, m_phi)
        
            dsig_dcos[m][i] = qabs[m][i]*A_v/(E_k[m] + m_e - kabs[m][i]*costheta[m][i])
        
            #if (i%100==0):
            #    print(m,i,A_v)
    
        func=inter.interp1d(costheta[m],dsig_dcos[m])

        def integrand(x):
            return func(x)
        I = quad(integrand, costheta[m][0], costheta[m][len(costheta[m])-1]) 

        Ne=64
        sig[m] = Ne*kappa**2*np.pi*alpha**2/(3*m_e*(E_k[m]**2-m_phi**2)**0.5)*I[0]
        sig[m] = sig[m]*(1.973*10**-11)**2 #MeV^-2 to cm^2
        sig[m] = sig[m]*10**24 #cm^2 to barn

    #sig = Calc_sig(dsig_dcos, sig, ppabs, pp, p, qabs, q, kabs, costheta, E_k, m_e, m_phi, kappa, alpha)

    # with open('/home/dfreitas/DPnum_Comp/Comp-xs2/CSxE_m=%f_k=%f.txt' %(m_phi, kappa), 'w') as archive:
    #     for i in range(0, len(E_k)):
    #         archive.write('%f   %.16f\n' %(E_k[i], sig[i]))

    plt.title('Cross-section')
    plt.plot(E_k, sig, c='blue', linewidth=2)
    plt.xlabel("$E_{k}[MeV]$")
    plt.ylabel("$\sigma~[barn/atom]$")
    plt.yscale('log')
    plt.xscale('linear')
    #plt.ylim(1E+2, 5e+9)
    #plt.xlim(1E-2,1)
    #plt.legend(loc = 'upper right')
    plt.grid(True)
    plt.show()
    #plt.savefig("/home/dfreitas/Results/CS_%dE_%dcos_k=%f_mphi=%f.png" %(len(E_k), len(q), kappa, m_phi))

main()