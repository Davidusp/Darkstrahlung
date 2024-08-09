import math
from math import isnan, nan
from os import kill
from tracemalloc import stop
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.interpolate as inter
import sys
import numba


# Calculate cross section considering rock (average of A=22 and Z=11)
@numba.jit(nopython=True, nogil=True)
def chi_calc(m_mu, E_mu, m_z, Z):
    A=22
    Z=11
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
def diff_cross_section(eps, alpha, m_z, E_mu, m_mu, Z, x):
    '''
        Calculate the differential cross section for DP production by muon Bremsstrahlung
        using the IWW approximation
    '''
    m_z = m_z/1000
    chi_IWW = chi_calc(m_mu, E_mu, m_z, Z)
    u_max = -m_z**2*(1-x)/x - m_mu**2*x
    u_min = -x*(E_mu*0.1)**2 + u_max #theta_max = 0.1 (E_mu(theta) = m_mu/theta)
    #print(u_max, u_min)
    dsig_dx_min = 2*eps**2*alpha**3*chi_IWW*np.sqrt((x**2 - (m_z/E_mu)**2))*(m_mu**2*x*(-2+2*x+x**2)-2*(3-3*x+x**2)*u_min)/(3*x*u_min**2)
    dsig_dx_max = 2*eps**2*alpha**3*chi_IWW*np.sqrt((x**2 - (m_z/E_mu)**2))*(m_mu**2*x*(-2+2*x+x**2)-2*(3-3*x+x**2)*u_max)/(3*x*u_max**2)
    dsig_dx = dsig_dx_max - dsig_dx_min

    return(dsig_dx)

#@numba.jit(nopython=True, nogil=True)
def ExtractData(text):
    '''
        This function separtes the data of an archive into diferent arrays
    '''    
    for i in range(len(text)):
        text[i] = text[i].split("\t")
        text[i] = [float(x) for x in text[i]]
    
    return(text)

#@numba.jit(nopython=True, nogil=True)
def Build_zenital(text, zenit, ratio, N):
    '''
        Build the arrays of zenital angles and the muon raios in this direction 
    '''
    for i in range(N):
        zenit[i]=np.arccos(text[i][0])
        ratio[i]=text[i][1]    

#@numba.jit(nopython=True, nogil=True)
def Build_general(text, list1, list2, N):

    for i in range(N):
        list1[i]=text[i][0]
        list2[i]=text[i][1]

def CalcOsc(Ei, theta, E, mDP, x, ratio, alpha, Nmuabove100, Na, epsilon, funcNaI, totaldecays):
    '''
    Calculate how many DPs are produced by muon Bremsstrahlung, and how many DPs would oscillate into photons inside COSINE LS. 
    Considered whole COSINE is made of NaI (need to correct this by considering crystals actual geometry). Also, considered photons 
    with energy below 2 MeV would be detected in the crystals with 100% efficiency. 
    '''

    N=40
    Z=11
    #m_z=np.logspace(-2.99,-1.301,30) #GeV
    #m_z = np.logspace(-5,-3.52,30)
    #E_mu = np.logspace(-1,3.7,2000) #GeV
    #print(np.diff(E_mu)[0])
    m_mu = 0.1057 #GeV/c^2
    alpha = 1/137
    L=200/(1.973*10**(-14)) #GeV^-1 # COSINE LS length
    me=0.511/1000 #GeV
    ne=9.433368*10**23*(1.973*10**(-14))**3 #For NaI (in GeV^3)
    mgamma=(4*np.pi*alpha*ne/me)**0.5

    dsig_dx=np.zeros((len(Ei),N))
    sig=np.zeros((len(epsilon),len(Ei),len(x)))
    ratio=np.zeros((len(Ei),N))

    f=open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Oscillation/ndecays_%s_oscillation_2MeV.txt' % "{:.2e}".format(mDP),'w')
    for k in range(len(epsilon)):
        for i in range(len(Ei)-1):
            if (i>0 and totaldecays[k]<1E-05): continue
            if (i>800): continue
            for t in range(len(theta)):
                for j in range(len(x)-1):
                    if ((700/np.cos(theta[t])-x[j]/2.65/100)>=0.0):
                        if (E[i][j]<1*10**2 or E[i][j]>5*10**3): continue
                        if (mDP*1e-3>=1E-6):
                            E_z = np.logspace(np.log10(mDP*1e-3), np.log10(0.002), N) #GeV
                        else:
                            E_z = np.logspace(np.log10(1E-6), np.log10(0.002), N) #GeV
                        for r in range(N):
                            ratio[i][r] = E_z[r]/(E[i][j]/1000)
                            if (r!=0):
                                dsig_dx[i][r] = diff_cross_section(epsilon[k], alpha, mDP, E[i][j]/1000, m_mu, Z, (ratio[i][r]+ratio[i][r-1])/2)
                                integr=(dsig_dx[i][r]+dsig_dx[i][r-1])*(ratio[i][r]-ratio[i][r-1])/2
                                sig[k][i][j]=integr*(1.973*10**(-16))**2*10**4*10**24*10**12 #pb

                                xsec=sig[k][i][j]*10**(-12)*10**(-24) #cm^2
                                Nbremss=Nmuabove100[i][t][j]*2.65*Na/22*np.diff(x)[0]/2.65*xsec*60*60*24*365.25
                        
                                Egamma=(E_z[r]+E_z[r-1])/2 #GeV
                                P=epsilon[k]**2*mDP**4*funcNaI(Egamma*1000)*L/((mDP**2-mgamma**2)**2+Egamma**2*funcNaI(Egamma*1000)**2)
                                totaldecays[k]=totaldecays[k]+P*Nbremss

            if (i%100==0 or i==800):
                print(k,i,totaldecays[k])
        print("{:.3e}".format(epsilon[k]),file=f,end=" ")
        print("{:.3e}".format(totaldecays[k]),file=f,end="\n")
    f.close()



def main():

    # Mass value given by user. It makes this script to be parallelizable.
    #mDP=float(sys.argv[1])
    mDP = 0.01  #MeV
    #print(mDP)

    # Import vertical muon flux from around 0.8 GeV to few TeV
    #arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonflux.tsv', 'r')
    arq = open('/home/dfreitas/DPnum_Comp/Muons/muonflux.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    p=np.zeros(N1)
    y=np.zeros(N1)
    E=np.zeros(N1)
    mu=0.105658 #GeV
    Na=6.02*10**23

    for i in range(N1):
        p[i]=texto1[i][0]
        y[i]=texto1[i][1]
        E[i]=(p[i]**2+mu**2)**0.5 # Since p>>mu, E~p
    
    # Interpolate muon flux in units of muons/cm^2/s/sr
    func=inter.interp1d(p,y/p**3)

    # Import muon stopping power in rock
    #arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muon_stoppower.tsv', 'r')
    arq = open('/home/dfreitas/DPnum_Comp/Muons/muon_stoppower.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    dEdx=np.zeros(N1) #Energy variation rate for the muon in each step
    Em=np.zeros(N1)   

    Build_general(texto1, Em, dEdx, N1)

    funcdEdx=inter.interp1d(Em,dEdx)

    # Import muon flux variation with the zenith angle for 1 GeV muons
    #arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonangle_1GeV.tsv', 'r')
    arq = open('/home/dfreitas/DPnum_Comp/Muons/muonangle_1GeV.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    zenite_1GeV=np.zeros(N1)
    ratio_1GeV=np.zeros(N1)

    Build_zenital(texto1, zenite_1GeV, ratio_1GeV, N1)

    func1GeV=inter.interp1d(zenite_1GeV,ratio_1GeV)

    # Import muon flux variation with the zenith angle for 10 GeV muons
    #arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonangle_10GeV.tsv', 'r')
    arq = open('/home/dfreitas/DPnum_Comp/Muons/muonangle_10GeV.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    zenite_10GeV=np.zeros(N1)
    ratio_10GeV=np.zeros(N1)

    Build_zenital(texto1, zenite_10GeV, ratio_10GeV, N1)

    func10GeV=inter.interp1d(zenite_10GeV,ratio_10GeV)

    # Import muon flux variation with the zenith angle for 100 GeV muons
    #arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonangle_100GeV.tsv', 'r')
    arq = open('/home/dfreitas/DPnum_Comp/Muons/muonangle_100GeV.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    zenite_100GeV=np.zeros(N1)
    ratio_100GeV=np.zeros(N1)

    Build_zenital(texto1, zenite_100GeV, ratio_100GeV, N1)

    func100GeV=inter.interp1d(zenite_100GeV,ratio_100GeV)

    # Import muon flux variation with the zenith angle for 1 TeV muons
    #arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonangle_1000GeV.tsv', 'r')
    arq = open('/home/dfreitas/DPnum_Comp/Muons/muonangle_1000GeV.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    zenite_1000GeV=np.zeros(N1)
    ratio_1000GeV=np.zeros(N1)

    Build_zenital(texto1, zenite_1000GeV, ratio_1000GeV, N1)

    func1000GeV=inter.interp1d(zenite_1000GeV,ratio_1000GeV)

        # Interpolate muon flux variation with zenith for different muon energies
    theta=np.linspace(3.0*np.pi/180,87*np.pi/180,15) # Zenith angles from 0 to 90 degrees, considering 3 degrees resolution
    #print(theta*180/np.pi)
    en=np.logspace(0,3,4)
    interen=np.zeros(4)
    Ei=np.logspace(3,6.69897,1000) #MeV # Muon energies from 1 GeV to few TeV
    #print((Ei[0]-np.diff(Ei)[0]/2)/1000)
    ratio=np.zeros((len(theta),len(Ei)))
    for t in range(len(theta)):
        interen[0]=func1GeV(theta[t])
        interen[1]=func10GeV(theta[t])
        interen[2]=func100GeV(theta[t])
        interen[3]=func1000GeV(theta[t])
        funcen=inter.interp1d(en,interen)
        for i in range(len(Ei)):
            if (Ei[i]<=1E+6):
                ratio[t][i]=funcen(Ei[i]/1000)
            else:
                ratio[t][i]=2.0*funcen(1E+6/1000) # Since the muon flux variation with zenith is only given to energies up to 1 TeV, consider that for energies above
                                                  # 1 TeV the muon flux is two times the muon flux for 1 TeV.

    x=np.linspace(0,(700/np.cos(theta[14]))/700*1855*100,5000) # Define step length in cm used in the muon attenuation and DP generated by Bremsstrahlung in the rock.

    # Calculate muon attenuation when travelling in the rock
    E=np.zeros((len(Ei),len(x)))
    for i in range(len(Ei)):
        
        for j in range(len(x)-1):
            if (j==0):
                E[i][j]=Ei[i]

            if (E[i][j]<=1.2): continue
            E[i][j+1]=E[i][j]-funcdEdx(E[i][j])*(x[j+1]-x[j])

    Acosine=(2*100)**2 #cm^2 COSINE-100 LS Area
    Nmu=np.zeros((len(theta),len(Ei)-1))

    # Check how many muons had a trajectory that would cross COSINE LS detector
    for i in range(len(Ei)-1):

        for t in range(len(theta)): 
            I = quad(func, (Ei[i]-np.diff(Ei)[i]/2)/1000, (Ei[i]+np.diff(Ei)[i]/2)/1000)
            sr=2*np.pi*(-np.cos(theta[t]+3.0*np.pi/180)+np.cos(theta[t]-3.0*np.pi/180))
            Nmu[t][i]=ratio[t][i]*I[0]*Acosine*sr #Muons/s

    # Check how many muons there are in each zenith angle, depth in rock, and with energy E[i][j]. 
    # It is required that the muon energy has to be greater than 100 MeV.
    Nmuabove100=np.zeros((len(Ei)-1,len(theta),len(x)))
    soma=0.0
    for i in range(len(Ei)-1):

        for t in range(len(theta)):
            j=0
            c=0
            while (E[i][j]>1.2 and np.cos(theta[t])*700*x[j]/1855/100<700):
                if (E[i][j]>1*10**2):
                    Nmuabove100[i][t][j]=Nmu[t][i]
    
                # Check how many muons hit COSINE. This value shoud be around 330 muons/day/cm^2.
                if(np.cos(theta[t])*700*x[j]/1855/100>695 and c==0):
                    soma=soma+Nmu[t][i]
                    c=c+1
                j=j+1

    alpha=1/137
    me=0.511/1000 #GeV

    epsilon=np.logspace(-10,-2,20)
    totaldecays=np.zeros(len(epsilon))

    # Import sodium mass attenuation
    arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Na_massatt.txt', 'r') #TROCAR MATERIAL!!!!!!!!!!!!!!!!!!
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    en_Na=np.zeros(N1)
    mu_Na=np.zeros(N1)

    for i in range(N1):
        en_Na[i]=texto1[i][0]
        mu_Na[i]=texto1[i][2]
    funcNa=inter.interp1d(en_Na,mu_Na)
   
    # Import iodine mass attenuation
    arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/I_massatt.txt', 'r') #TROCAR MATERIAL!!!!!!!!!!!!!!!!!!
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    en_I=np.zeros(N1)
    mu_I=np.zeros(N1)

    for i in range(N1):
        en_I[i]=texto1[i][0]
        mu_I[i]=texto1[i][2]
    funcI=inter.interp1d(en_I,mu_I)

    en=[]
    enfinal=[]
    for i in range(len(en_Na)):
        en.append(en_Na[i])
    for i in range(len(en_I)):
        en.append(en_I[i])
    for i in range(len(en)):
        for j in range(i+1,len(en)):
            if (en[i]==en[j]):
                en[i]=0.0
    for i in range(len(en)):
        if (en[i]!=0.0):
            enfinal.append(en[i])  

    att_NaI=np.zeros(len(enfinal))
    for i in range(len(enfinal)):
        att_NaI[i]=(funcNa(enfinal[i])+funcI(enfinal[i]))/2
    funcNaI=inter.interp1d(enfinal,att_NaI*3.67*(1.973*10**(-14))) #GeV, with energy given in MeV # Mass attenuatioin for NaI.


    CalcOsc(Ei, theta, E, mDP, x, ratio, alpha, Nmuabove100, Na, epsilon, funcNaI, totaldecays)    

main()