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
def chi_calc(m_mu, E_mu, m_z):
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
def diff_cross_section(eps, alpha, m_z, E_mu, m_mu, x):
    '''
        Calculate the differential cross section for DP production by muon Bremsstrahlung
        using the IWW approximation
    '''
    m_z = m_z/1000
    chi_IWW = chi_calc(m_mu, E_mu, m_z)
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

def ExtractData1(text):
    '''
        This function separtes the data of an archive into diferent arrays
    '''    
    for i in range(len(text)):
        text[i] = text[i].split("  ")
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

def InterpMassAtt(element):

    # Import sodium mass attenuation
    arq = open('/home/dfreitas/DPnum_Comp/Oscillation/Mass_att(Exmu)_%s.txt' %(element), 'r') 
    #arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Oscillation/Mass_att(Exmu)_%s.txt' %(element), 'r') 
    texto = arq.readlines()
    arq.close()

    texto = ExtractData1(texto)

    N1=len(texto)
    en=np.zeros(N1)
    mu=np.zeros(N1)

    for i in range(N1):
        en[i]=texto[i][0]
        mu[i]=texto[i][2]
    func=inter.interp1d(en,mu)

    return(func, en)

#@numba.jit(nopython=True, nogil=True)
def GatherEnergies(en_A,en_B):

    en=[]
    enfinal=[]
    for i in range(len(en_A)):
        en.append(en_A[i])
    for i in range(len(en_B)):
        en.append(en_B[i])
    for i in range(len(en)):
        for j in range(i+1,len(en)):
            if (en[i]==en[j]):
                en[i]=0.0
    for i in range(len(en)):
        if (en[i]!=0.0):
            enfinal.append(en[i])  
    
    return(enfinal)

def CalcOsc(Ei, theta, E, mDP, x, alpha, Nmuabove100, Na, epsilon, funcrock):
    '''
    Calculate how many DPs are produced by muon Bremsstrahlung, and how many DPs would oscillate into photons inside the rock above. 
    COSINE. Considered whole COSINE is made of NaI (need to correct this by considering crystals actual geometry). Also, considered 
    photons with energy below 2 MeV would be detected in the crystals with 100% efficiency. 
    '''

    N=50
    Z=11
    
    m_mu = 0.1057 #GeV/c^2
    me=0.511/1000 # electron mass (GeV)
    ne=2.65/2/(1.67262192*10**(-24))*(1.973*10**(-14))**3 #Electron number density for rock (GeV^3)
    mgamma=(4*np.pi*alpha*ne/me)**0.5  #effective photon mass (GeV)
    L = np.diff(x)[0]/2.65 #step length inside the rock (cm)

    dsig_dx=np.zeros((len(Ei),N))
    sig=np.zeros((len(Ei),len(x)))
    ratio=np.zeros((len(Ei),N))

    NDPs = np.zeros((len(theta),len(x),N))
    NOsc = np.zeros((len(theta),len(x),N))

    if (mDP*1e-3>=1E-6):
        E_z = np.linspace(mDP*1e-3, 0.1, N) #GeV
    else:
        E_z = np.linspace(1E-6, 0.1, N) #GeV
  
    #for k in range(len(epsilon)):
    for i in range(len(Ei)-1):
        if i%10 == 0:
            print(i)
        #if (i>0 and totaldecays<1E-05): continue
        if (i>=900): continue
        for t in range(len(theta)):
            for j in range(len(x)-1):
               
                if ((700/np.cos(theta[t])-x[j]/2.65/100)>=0.0):
                    if (E[i][j]<1*10**2): continue
                    
                    for r in range(N):
                        ratio[i][r] = E_z[r]/(E[i][j]/1000)
                        if (r!=0):
                            dsig_dx[i][r] = diff_cross_section(epsilon, alpha, mDP, E[i][j]/1000, m_mu, (ratio[i][r]+ratio[i][r-1])/2)
                            integr=(dsig_dx[i][r]+dsig_dx[i][r-1])*(ratio[i][r]-ratio[i][r-1])/2
                            sig[i][j]=integr*(1.973*10**(-16))**2*10**4*10**24*10**12 #pb
                            xsec=sig[i][j]*10**(-12)*10**(-24) #cm^2
                            Nbremss=Nmuabove100[i][t][j]*2.65*Na/22*L*xsec*60*60*24*365.25
                            NDPs[t][j][r] += Nbremss
                    
                            # Egamma=(E_z[r]+E_z[r-1])/2 #GeV
                            # P = epsilon**2*(mDP*1e-03)**4*funcrock(Egamma*1000)*L*5.06e13/(((mDP*1e-03)**2-mgamma**2)**2+Egamma**2*funcrock(Egamma*1000)**2)

                            # if j>0:
                            #     DPshere = (np.sum(NDPs[t][:j+1,r]) - np.sum(NOsc[t][:j,r]))
                            #     NOsc[t][j][r] += round(DPshere,3)*P

                            #     if np.sum(NDPs[t][:j+1,r]) < DPshere or P>1:
                            #         print(P, Nbremss, DPshere, i, t, j, funcrock(Egamma*1000), L, Egamma, mgamma)
                            #         #exit()
                            # else:
                            #     NOsc[t][j][r] += P*NDPs[t][j][r]

                            #     if P>1:
                            #         print(P, Nbremss, DPshere, i, t, j, funcrock(Egamma*1000), L, Egamma, mgamma)
                            #         #exit()
                            
    for t in range(len(theta)):
        for j in range(len(x)-1):
               
            if ((700/np.cos(theta[t])-x[j]/2.65/100)>=0.0):
                if (E[i][j]<1*10**2): continue
                    
                for r in range(N):
                    if (r!=0):
                        Egamma=(E_z[r]+E_z[r-1])/2 #GeV
                        P = epsilon**2*(mDP*1e-03)**4*funcrock(Egamma*1000)*L*5.06e13/(((mDP*1e-03)**2-mgamma**2)**2+Egamma**2*funcrock(Egamma*1000)**2)

                        if j>0:
                            
                            DPshere = (np.sum(NDPs[t][:j+1,r]) - np.sum(NOsc[t][:j,r]))
                            NOsc[t][j][r] += P*round(DPshere,3)

                            if np.sum(NDPs[t][:j+1,r]) < DPshere or P>1:
                                print(P, Nbremss, DPshere, i, t, j, funcrock(Egamma*1000), L, Egamma, mgamma)
                                exit()
                        else:
                            NOsc[t][j][r] += P*NDPs[t][j][r]                

    return(NOsc)

def main():

    # Mass value given by user. It makes this script to be parallelizable.
    mDP=float(sys.argv[1])
    epsilon=float(sys.argv[2])
    # mDP = 0.000121 #MeV
    # epsilon = 0.00133352143
    # #print(mDP)

    # Import vertical muon flux from around 0.8 GeV to few TeV
    #arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonflux.tsv', 'r')
    arq = open('/home/dfreitas/DPnum_Comp/Muons/muonflux.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    p=np.zeros(N1)
    y=np.zeros(N1)
    #E=np.zeros(N1)
    mu=0.105658 #GeV
    Na=6.02*10**23

    for i in range(N1):
        p[i]=texto1[i][0]
        y[i]=texto1[i][1]
        #E[i]=(p[i]**2+mu**2)**0.5 # Since p>>mu, E~p
    
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

    for i in range(len(Ei)-1):

        for t in range(len(theta)):
            j=0
            #c=0
            while (E[i][j]>1.2 and np.cos(theta[t])*700*x[j]/1855/100<700):
                if (E[i][j]>1*10**2):
                    Nmuabove100[i][t][j]=Nmu[t][i]
    
                # Check how many muons hit COSINE. This value shoud be around 330 muons/day/cm^2.
                # if(np.cos(theta[t])*700*x[j]/1855/100>695 and c==0):
                #     soma=soma+Nmu[t][i]
                #     c=c+1
                j=j+1

    alpha=1/137

    (funcNa, en_Na)=InterpMassAtt("Na")
    (funcSi, en_Si)=InterpMassAtt("Si")
    (funcAl, en_Al)=InterpMassAtt("Al")
    (funcO, en_O)=InterpMassAtt("O")
    (funcFe, en_Fe)=InterpMassAtt("Fe")
    (funcCa, en_Ca)=InterpMassAtt("Ca")
    (funcK, en_K)=InterpMassAtt("K")
    (funcMg, en_Mg)=InterpMassAtt("Mg")
    (funcTi, en_Ti)=InterpMassAtt("Ti")
    
    en1 = GatherEnergies(en_Na, en_Si)
    en2 = GatherEnergies(en1, en_Al)  
    en3 = GatherEnergies(en2, en_O)
    en4 = GatherEnergies(en3, en_Fe)
    en5 = GatherEnergies(en4, en_K)
    en6 = GatherEnergies(en5,en_Ca)
    en7 = GatherEnergies(en6,en_Ti)
    enfinal = GatherEnergies(en7, en_Mg)

    att_rock=np.zeros(len(enfinal))
    for i in range(len(enfinal)): #http://www.physicalgeography.net/fundamentals/10d.html
        att_rock[i]= 0.0283*funcNa(enfinal[i])+0.2772*funcSi(enfinal[i])+0.466*funcO(enfinal[i])+0.0813*funcAl(enfinal[i])
        +0.05*funcFe(enfinal[i])+0.0363*funcCa(enfinal[i])+0.0259*funcK(enfinal[i])+0.0209*funcMg(enfinal[i])+0.0148*funcTi(enfinal[i])
    funcrock=inter.interp1d(enfinal, att_rock*2.65*(1.973*10**(-14)), kind='linear', fill_value="extrapolate") #GeV, with energy given in MeV # Mass attenuatioin for the rock.

    NOsc = CalcOsc(Ei, theta, E, mDP, x, alpha, Nmuabove100, Na, epsilon, funcrock)    

    with open('/home/dfreitas/DPnum_Comp/Oscillation/Num_right/noscRock_m=%.13f_k=%.11f.txt' % (mDP,epsilon), 'w') as archive:
        archive.write('%s   %s\n' %("{:e}".format(epsilon),"{:.20e}".format(np.sum(NOsc))))

    # with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Oscillation/noscRock_m=%f_k=%f.txt' % (mDP,epsilon), 'w') as archive:
    #     archive.write('%s   %s\n' %("{:e}".format(epsilon),"{:.20e}".format(np.sum(NOsc))))

main()
