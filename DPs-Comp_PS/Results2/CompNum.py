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

def CalcBremss(Ei, theta, E, mDP, N, x, kappas, ratio, dsig_dx, alpha, sig, Nmuabove100, m_mu, Z, Na):
    '''
    Function to compute the number NDPs_tot of muon Bremsstralung at the rock
    producing DPs at each step j and angle t, for each coupling parameter and a
    DP mass mDP.
    '''

    #f=open('/home/davidc/Documents/Master\'s Analisys/Parameter space/nbremss_%s_%s.txt' %("{:.3e}".format(mDP),"{:.3e}".format(epsilon)),'w')
    
    NDPs_tot = np.zeros((len(theta), len(x))) #number of DPs generated in each angle and step
    #print(f)
    #totalbremss = 0

    for i in range(len(Ei)-1):
    
        if (i>=900): continue
        for t in range(0,len(theta)):
            
            #cont = 0
            for j in range(0,len(x)-1):

                if ((700/np.cos(theta[t])-x[j]/2.65/100)>=0.0):
                    if (E[i][j]<1*10**2): continue
                    if (mDP>1E-06):
                        E_z = np.logspace(np.log10(mDP*1e-3), np.log10(0.1), N) #GeV
                    
                    else:
                        E_z = np.logspace(np.log10(1E-06), np.log10(0.1), N) #GeV
                    for r in range(N):
                        ratio[i][r] = E_z[r]/(E[i][j]/1000)
            
                        if (r!=0):
                            dsig_dx[i][r] = diff_cross_section(kappas, alpha, mDP, E[i][j]/1000, m_mu, Z, (ratio[i][r]+ratio[i][r-1])/2)
                            integr=(dsig_dx[i][r]+dsig_dx[i][r-1])*(ratio[i][r]-ratio[i][r-1])/2
                            sig[i][j]=integr*(1.973*10**(-16))**2*10**4*10**24*10**12 #pb
                            xsec=sig[i][j]*10**(-12)*10**(-24) #cm^2
                            Nbremss=Nmuabove100[i][t][j]*2.65*Na/22*np.diff(x)[0]/2.65*xsec*60*60*24*365.25 # num of DP's per year
                            NDPs_tot[t][j] += Nbremss
                            #totalbremss += Nbremss

            #print(i,totalbremss)        
    #print("{:.8e}".format(totalbremss),file=f,end="\n")

    #f.close()
    return(NDPs_tot)

def CalcComp(theta, mDP, x, epsilon, Na, kappas, NDPs_tot):
    '''
    Function to estimate the number of Compton-like scatterings occuring for each angle
    and step inside the rock, and for each coupling parameter and DP mass 
    '''

    NComps_res = np.zeros((len(theta), len(x))) #Number of Compton-like scatterings in each angle and step

    # with open('/home/davidc/Documents/Master\'s Analisys/Results Compton/CSxE_m=%s_k=%s00.txt' %(mDP, epsilon), 'r') as archive:
    #         arq = archive.read() 

    with open('/home/dfreitas/DPnum_Comp/Comp-xs2/CSxE_m=%f_k=%s00.txt' %(mDP, epsilon), 'r') as archive:
            arq = archive.read()

    lines = arq.split("\n")
    lines = lines[:-1]

    for row in range(0, len(lines)):
        lines[row] = lines[row].split("   ")
        lines[row] = [float(x) for x in lines[row]]
        #print(lines)
    #arq = arq.split("\t")
    #print(lines)

    sigComp = np.zeros(len(lines))
    E_k = np.zeros(len(lines))

    for l in range(0, len(lines)):
        E_k[l] = lines[l][0]
        sigComp[l] = lines[l][1]

    CScompfunc = inter.interp1d(E_k, sigComp, kind = 'cubic', fill_value="extrapolate")

    E_k = np.logspace(mDP+0.01*mDP,np.log10(100), len(x))

    sigComp = (kappas/epsilon)**2*CScompfunc(E_k)
    # print(sigComp)

    for t in range(0,len(theta)):
            
        for j in range(0,len(x)-1):

            if ((700/np.cos(theta[t])-x[j]/2.65/100)>=0.0):

                if j>0:
                    NComps = (np.sum(NDPs_tot[t][:j+1]) - np.sum(NComps_res[t][:j]))*2.65*Na/22*np.diff(x)[0]/2.65*sigComp[len(sigComp)-1-j]*1e-24 # num of Comptons per year
                    #print(np.sum(NDPs_tot[t][:j+1]), np.sum(NComps_res[k][t][:j]))
                    NComps_res[t][j] += NComps

                else:
                    NComps = NDPs_tot[t][j]*2.65*Na/22*np.diff(x)[0]/2.65*sigComp[len(sigComp)-1-j]*1e-24 # num of Comptons per year
                    NComps_res[t][j] += NComps

    return(NComps_res)

def main():
    # Mass and epsilon values given by user. It makes this script to be parallelizable.
    mDP=float(sys.argv[1]) #MeV
    kappas=float(sys.argv[2])
    
    #mDP = 0.323746 #MeV
    epsilon = 1e-4
    #kappas = 1e-10

    #print(mDP,epsilon)

    # Import vertical muon flux from around 0.8 GeV to few TeV
    #arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonflux.tsv', 'r')
    arq = open('/home/dfreitas/DPnum_Comp/Muons/muonflux.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    p=np.zeros(N1)    #muon momentum
    y=np.zeros(N1)
    #E=np.zeros(N1)   #muon energy
    mu=0.105658 #GeV  #muon mass
    Na=6.02*10**23    #Avogadro number

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

    # Check how many muons there are in each zenith angle, depth in rock, and with energy E[i][j]. It is required that the muon energy has to be greater than 100 MeV.
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

    #print(soma/4*3600*24)

    alpha=1/137
    me=0.511/1000 #GeV
    m_mu=105.7/1000 #GeV

    N=50
    Z=11
    dsig_dx=np.zeros((len(Ei),N))
    sig=np.zeros((len(Ei),len(x)))
    ratio=np.zeros((len(Ei),N))

    # Calculate how many DPs are produced by muon Bremsstrahlung.

    NDPs_tot = CalcBremss(Ei, theta, E, mDP, N, x, kappas, ratio, dsig_dx, alpha, sig, Nmuabove100, m_mu, Z, Na)

    NComps_res = CalcComp(theta, mDP, x, epsilon, Na, kappas, NDPs_tot)

    # with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/DPtot_Comp_m=%f_k=%.11f.txt' %(mDP,kappas), 'w') as archive:
    #     totalDPs = NDPs_tot - NComps_res
    #     archive.write('%.11f   %.16f\n' %(kappas, np.sum(totalDPs)))

    with open('/home/dfreitas/DPnum_Comp/totalDPs/Results2/DPtot_Comp_m=%.8f_k=%.11f.txt' %(mDP,kappas), 'w') as archive:
        totalDPs = NDPs_tot - NComps_res
        archive.write('%.11f   %.16f\n' %(kappas, np.sum(totalDPs)))

    #print("total DPs = \n", totalDPs)
    #print("Number of Bremsstrahlung producing DPs: \n", np.sum(NDPs_tot))
    #print("\nNumber of Compton-like scatterings: \n", "{:e}".format(np.sum(NComps_res)))
    #print("\nNumero total de Comptons: %d." %(np.sum(totalDPs)))
    #print("\nPorcentagem de DP's que sofrem Compton: %d" %(np.sum(NComps_res)/np.sum(NDPs_tot)))

main()