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

# Calculate cross section considering rock (average of A=22 and Z=11)
@numba.jit(nopython=True, nogil=True)
def chi_calc(m_mu, E_mu, m_z, Z):
    A=22
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

def ExtractData1(text):
    '''
        This function separtes the data of an archive into diferent arrays
    '''    
    for i in range(len(text)):
        text[i] = text[i].split("  ")
        text[i] = [float(x) for x in text[i]]
    
    return(text)

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

def CalcEvents(Ei, theta, E, mDP, x, kappas, Nmuabove100, Na, epsilon, funcrock, funcNaI):
    '''
    Function to compute the number NDPs_tot of muon Bremsstralung at the rock
    producing DPs at each step j and angle t, for each coupling parameter and a
    DP mass mDP.
    '''
    alpha1=1/137
    me=0.511/1000 #electron mass (GeV)
    m_mu=105.7/1000 #muon mass (GeV)
    L = np.diff(x)[0]/2.65 #step length inside the rock (cm)

    N=100
    Z=11
    dsig_dx=np.zeros((len(Ei),N))
    sig=np.zeros((len(Ei),len(x)))
    ratio=np.zeros((len(Ei),N))

    #FIRST WE COMPUTE THE PROCESSES IN THE ROCK
    NDPs_tot = np.zeros((len(theta), len(x), N)) #Number of DPs generated in each angle, step and energy
    NComps_res = np.zeros((len(theta), len(x), N)) #Number of Comptons in each angle, step and energy
    NAbs_res = np.zeros((len(theta), len(x), N)) #Number of absorption processes in each angle, step and energy
    NOsc_res = np.zeros((len(theta), len(x), N)) #Number of oscillations in each angle, step and energy

    q = 0.30282212088 #Elementary charge in natural units
    alpha = q**2/(4*np.pi) #Fine structure constant
    alphaline = (q*kappas)**2/(4*np.pi)

    ne=2.65/2/(1.67262192*10**(-24))*(1.973*10**(-14))**3 #Electron number density for rock (GeV^3)
    mgamma=(4*np.pi*alpha1*ne/me)**0.5  #effective photon mass (GeV)

    #Compton-like scaterring Cross-Section
    with open('/home/dfreitas/DPnum_Comp/Comp-xs2/CSxE_m=%f_k=%s00.txt' %(mDP, epsilon), 'r') as archive1:
            arqComp = archive1.read()

    # with open('/home/davidc/Documents/Master\'s Analisys/Results Compton/Cross Section 2/CSxE_m=%f_k=%s00.txt' %(mDP, epsilon), 'r') as archive1:
    #         arqComp = archive1.read() 
        
    linesComp = arqComp.split("\n")
    linesComp = linesComp[:-1]

    for row in range(0, len(linesComp)):
        linesComp[row] = linesComp[row].split("   ")
        linesComp[row] = [float(x) for x in linesComp[row]]
    
    sigComparq = np.zeros(len(linesComp))
    E_kComp = np.zeros(len(linesComp))

    for l in range(0, len(linesComp)):
        E_kComp[l] = linesComp[l][0]
        sigComparq[l] = (Z/64)*linesComp[l][1]*(1e-24)*(5.06e13)**2 #barn/atom->cm²/atom->GeV⁻²/atom
    
    CScompfunc = inter.interp1d(E_kComp, sigComparq, kind = 'cubic', fill_value="extrapolate")

    #Absorption Cross-Section
    # with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/PEXSRock.txt', 'r') as archive2:
    #         arqAbs = archive2.read()
    with open('/home/dfreitas/DPnum_Comp/Absorption/PEXSRock.txt', 'r') as archive2:
            arqAbs = archive2.read()  

    linesAbs = arqAbs.split("\n")
    linesAbs = linesAbs[:-1]
    
    for row in range(0, len(linesAbs)):
        linesAbs[row] = linesAbs[row].split(" ")
        linesAbs[row] = [float(x) for x in linesAbs[row]]

    sigAbs = np.zeros(len(linesAbs))
    E_kAbs = np.zeros(len(linesAbs))

    for l in range(0, len(linesAbs)):
        E_kAbs[l] = linesAbs[l][0]
        if mDP < E_kAbs[l]:
            sigAbs[l] = (linesAbs[l][1]*(5.06e13)**2/np.sqrt(1-mDP**2/E_kAbs[l]**2))*(alphaline/alpha)*(22/Na) #GeV⁻²/atom (1/(Na/A)))

    CSabsfunc = inter.interp1d(E_kAbs, sigAbs, kind = 'linear', fill_value="extrapolate")

    if (mDP*1e-3>1E-06):
        E_z = np.linspace(mDP*1e-3, 0.0002, N) # DP energy (GeV)
                    
    else:
        E_z = np.linspace(1E-06, 0.0002, N) #DP energy (GeV)

    for i in range(len(Ei)-1):
        
        if (i>=900): continue
        for t in range(0,len(theta)):
            
            for j in range(0,len(x)-1):

                if ((700/np.cos(theta[t])-x[j]/2.65/100)>=0.0):
                    if (E[i][j]<1*10**2): continue
                    for r in range(N):
                        ratio[i][r] = E_z[r]/(E[i][j]/1000)
            
                        if (r!=0):

                            #Bremsstrahlung
                            dsig_dx[i][r] = diff_cross_section(kappas, alpha1, mDP, E[i][j]/1000, m_mu, Z, (ratio[i][r]+ratio[i][r-1])/2)
                            integr=(dsig_dx[i][r]+dsig_dx[i][r-1])*(ratio[i][r]-ratio[i][r-1])/2
                            sig[i][j]=integr*(1.973*10**(-16))**2*10**4*10**24*10**12 #pb
                            xsec=sig[i][j]*10**(-12)*10**(-24) #cm^2
                            Nbremss=Nmuabove100[i][t][j]*2.65*Na/22*L*xsec*60*60*24*365.25 # num of DP's per year
                            NDPs_tot[t][j][r] += Nbremss
              
    for t in range(0,len(theta)): #angles
        for j in range(0,len(x)-1): #steps inside rock
            if ((700/np.cos(theta[t])-x[j]/2.65/100)>=0.0):
                    
                sigComp = (kappas/epsilon)**2*CScompfunc(1e03*E_z) #Compton cross-section
                sigAbs = CSabsfunc(1e03*E_z)  #Absorption cross section

                for r in range(0, N):

                    if r!= 0:

                        Egamma=(E_z[r]+E_z[r-1])/2 #GeV
                        P=kappas**2*(mDP*1e-03)**4*funcrock(Egamma*1000)*L*5.06e13/(((mDP*1e-03)**2-mgamma**2)**2+Egamma**2*funcrock(Egamma*1000)**2)
                        if P<=1e-15:
                            P=0
                        if P>1:
                            P=1

                        if j>0:
                            DPshere = (np.sum(NDPs_tot[t][:j+1,r]) - np.sum(NComps_res[t][:j,r]) 
                                        - np.sum(NAbs_res[t][:j,r]) - np.sum(NOsc_res[t][:j,r]))
                            
                            if np.abs(DPshere) <= 1e-10:
                                DPshere = 0

                            #Compton
                            lambComp = 1/(ne*sigComp[r])
                            #NComps = (round(DPshere, 4)-NOsc)*(1 - np.exp(-L*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                            NComps = (DPshere)*(1 - np.exp(-L*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                            if np.abs(NComps) <=1e-10:
                                NComps = 0
                            NComps_res[t][j][r] += NComps

                            if np.sum(NDPs_tot[t][:j+1,r]) < np.sum(NComps_res[t][:j,r]):
                                print("j>0, Comp", NComps, lambComp, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NComps_res[t][:j,r]), np.exp(-L*5.06e13/(lambComp)), NDPs_tot[t][j][r], NComps_res[t][j][r], sigComp[r],t,j)
                                #exit()
                                NComps_res[t][j][r] -= np.abs(np.sum(NComps_res[t][:j,r]) - np.sum(NDPs_tot[t][:j+1,r]))
                            
                            #Absorption
                            lambAbs = 1/(ne*sigAbs[r])
                            #NAbs = (round(DPshere, 4)-NOsc-NComps)*(1 - np.exp(-L*5.06e13/(lambAbs)))  #num of Absorptions/year for this DP energy
                            NAbs = (DPshere-NComps)*(1 - np.exp(-L*5.06e13/(lambAbs)))  #num of Absorptions/year for this DP energy
                            if np.abs(NAbs) <=1e-10:
                                NAbs = 0
                            NAbs_res[t][j][r] += NAbs
                            
                            if np.sum(NDPs_tot[t][:j+1,r]) < np.sum(NAbs_res[t][:j,r]):
                                print("j>0, Abs", NAbs, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NAbs_res[t][:j,r]), np.exp(-L*5.06e13/(lambAbs)), sigAbs[r],t,j,E_z[r])
                                #exit()
                                NAbs_res[t][j][r] -= np.abs(np.sum(NAbs_res[t][:j,r])-np.sum(NDPs_tot[t][:j+1,r]))

                            #Oscillation
                            #NOsc = round(DPshere,3)*P
                            NOsc = (DPshere-NComps-NAbs)*P
                            if np.abs(NOsc) <=1e-10:
                                NOsc = 0
                            NOsc_res[t][j][r] += NOsc 

                            if np.sum(NDPs_tot[t][:j+1,r]) < np.sum(NOsc_res[t][:j,r]):
                                print("j>0, Osc", P, NOsc, DPshere, j)
                                #exit()
                                NOsc_res[t][j][r] -= np.abs(np.sum(NOsc_res[t][:j,r])-np.sum(NDPs_tot[t][:j+1,r]))
                            
                            if DPshere<0:
                                print("j>0, Pcomp = %f, Pabs = %f, Posc = %f, DPshere = %f, NComps = %f, NAbs = %f, NOsc = %f, NDPsSum = %f," 
                                  "NCompsSum = %f, NAbsSum = %f,NOscSum = %f, j=%d, r=%d" %(1-np.exp(-L*5.06e13/(lambComp)), 
                                   1-np.exp(-L*5.06e13/(lambAbs)), P, DPshere, NComps, NAbs, NOsc, np.sum(NDPs_tot[t][:j+1,r]), 
                                   np.sum(NComps_res[t][:j,r]), np.sum(NAbs_res[t][:j,r]),np.sum(NOsc_res[t][:j,r]) , j, r))
 
                                exit()

                        else:

                            #Compton
                            lambComp = 1/(ne*sigComp[r])
                            NComps = (NDPs_tot[t][j][r])*(1 - np.exp(-L*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                            if np.abs(NComps) <=1e-10:
                                NComps = 0
                            NComps_res[t][j][r] += NComps

                            if NComps_res[t][j][r] > NDPs_tot[t][j][r] or NComps > NDPs_tot[t][j][r]:
                                print("j=0, Comp",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambComp)), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambComp)), NDPs_tot[t][j][r])
                                #exit()
                                NComps_res[t][j][r] -= np.abs(NDPs_tot[t][j][r]-NComps_res[t][j][r])

                            #Absorption
                            lambAbs = 1/(ne*sigAbs[r])
                            NAbs = (NDPs_tot[t][j][r]-NComps)*(1 - np.exp(-L*5.06e13/(lambAbs))) # num of absorptions per year for this DP energy
                            if np.abs(NAbs) <=1e-10:
                                NAbs = 0
                            NAbs_res[t][j][r] += NAbs

                            if NAbs_res[t][j][r] > NDPs_tot[t][j][r] or NAbs > NDPs_tot[t][j][r]:
                                print("j=0, Abs",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambAbs)), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambAbs)), NDPs_tot[t][j][r])
                                #exit()
                                NAbs_res[t][j][r] -= np.abs(NDPs_tot[t][j][r]-NAbs_res[t][j][r])

                            #Oscillation
                            NOsc = (NDPs_tot[t][j][r]-NAbs-NComps)*P
                            if np.abs(NOsc) <=1e-10:
                                NOsc = 0
                            NOsc_res[t][j][r] += NOsc

                            if NDPs_tot[t][j][r] < NOsc:
                                print("j=0, Osc", P, NOsc, DPshere2)
                                #exit()
                                NOsc_res[t][j][r] -= np.abs(NOsc - NDPs_tot[t][j][r])

                            if NDPs_tot[t][j][r] < (NComps_res[t][j][r]+NAbs+NOsc):
                                print("j=0", np.exp(-L*5.06e13/(lambComp)), np.exp(-L*5.06e13/(lambAbs)), P, NDPs_tot[t][j][r], NComps, NAbs, NOsc)

    NDPs_res = np.zeros((len(theta),N))
    #print(NDPs_tot)
    for t in range(0,len(theta)): #angles
        for r in range(0, N):
            if r!= 0:
                NDPs_res[t][r] = (np.sum(NDPs_tot[t][:,r])-(np.sum(NComps_res[t][:,r])+np.sum(NAbs_res[t][:,r])+np.sum(NOsc_res[t][:,r])))

    #print(NDPs_res)

    # for t in range(0,len(theta)): #angles
    #     for j in range(0,len(x)-1): #steps inside rock
    #         if ((700/np.cos(theta[t])-x[j]/2.65/100)>=0.0):
    #             for r in range(0, N):
    #                 if r!= 0:
    #                     NDPs_tot[t][j][r] = round(NDPs_tot[t][j][r], 1)

    #NOW WE COMPUTE THE INTERACTIONS WITH THE CRYSTAL
    xdetect = np.linspace(0,(0.4/np.cos(theta[14]))*100,800) #positions inside the crystals (cm)
    NComps_tot = np.zeros((len(theta), len(xdetect), N)) #Number of Comptons in each angle, step and energy
    NAbs_tot = np.zeros((len(theta), len(xdetect), N)) #Number of absorptions in each angle, step and energy
    NOsc_tot = np.zeros((len(theta), len(xdetect), N)) #Number of oscillations in each angle, step and energy
    rho_NaI = 3.67  #NaI density (g/cm³)
    Ldetec = np.diff(xdetect)[0] #(cm)
    ne_NaI=rho_NaI/(149.89/Na)*(1.973*10**(-14))**3 #Electron number density for NaI (GeV³)
    mgamma_NaI=(4*np.pi*alpha1*ne_NaI/me)**0.5    #effective photon mass

    CScompfuncNaI = inter.interp1d(E_kComp, sigComparq, kind = 'cubic', fill_value="extrapolate")

    # with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/PeXSNaI.txt', 'r') as archive3:
    #         arqAbs1 = archive3.read() 
    with open('/home/dfreitas/DPnum_Comp/Absorption/PeXSNaI.txt', 'r') as archive3:
            arqAbs1 = archive3.read()

    linesAbs1 = arqAbs1.split("\n")
    linesAbs1 = linesAbs1[:-1]
    
    for row in range(0, len(linesAbs1)):
        linesAbs1[row] = linesAbs1[row].split(" ")
        linesAbs1[row] = [float(x) for x in linesAbs1[row]]

    sigAbsNaI = np.zeros(len(linesAbs1))
    E_kAbs = np.zeros(len(linesAbs1))

    for l in range(0, len(linesAbs1)):
        E_kAbs[l] = linesAbs1[l][0]
        if mDP < E_kAbs[l]:
            sigAbsNaI[l] = (linesAbs1[l][1]*(5.06e13)**2/np.sqrt(1-mDP**2/E_kAbs[l]**2))*(alphaline/alpha)*(149.89/Na) #GeV⁻²/atom (1/(Na/A))) #barns/atom
    
    CSabsfuncNaI = inter.interp1d(E_kAbs, sigAbsNaI, kind = 'linear', fill_value="extrapolate")

    if (mDP>1E-03):
        E_z = np.linspace(mDP, 0.2, N) #MeV
                    
    else:
        E_z = np.linspace(1E-03, 0.2, N) #MeV

    for t in range(0,len(theta)):
        for j in range(0,len(xdetect)-1):
            if ((0.4/np.cos(theta[t])-xdetect[j]/100)>=0.0):
                    
                sigCompNaI = (kappas/epsilon)**2*CScompfuncNaI(E_z) #Compton cross-section
                sigAbsNaI = CSabsfuncNaI(E_z)  #Absorption cross section

                for r in range(N):

                    if j>0:
                        DPshere2 = NDPs_res[t][r]
                        
                        if np.abs(DPshere2) <= 1e-12:
                            DPshere2 = 0
                        
                        #Compton
                        lambComp = 1/(ne_NaI*sigCompNaI[r])
                        #NComps = (round(DPshere2,2))*(1 - np.exp(-Ldetec*5.06e13/(lambComp))) # num of Compton/year inside NaI crystals
                        NComps = (DPshere2)*(1 - np.exp(-Ldetec*5.06e13/(lambComp))) # num of Compton/year inside NaI crystals
                        if np.abs(NComps) <=1e-12:
                            NComps = 0
                            #print("sigComp = %f, P = %f, E = %f" %(sigCompNaI[r], 1 - np.exp(-Ldetec*5.06e13/(lambComp)), E_z[r]))
                        # else:
                        #     print("NComps = %.16f" %(NComps))
                            #exit()

                        NComps_tot[t][j][r] += NComps

                        if np.sum(NDPs_tot[t][:j+1,r]) < np.sum(NComps_tot[t][:j,r]):
                            print("j>0, CompNaI", NComps, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NComps_tot[t][:j,r]), DPshere2, sigCompNaI[r],t,j,E_z[r])
                            exit()
                            #NComps_tot[t][j][r] -= np.abs(np.sum(NComps_tot[t][:j,r])-np.sum(NDPs_tot[t][:j+1,r]))

                        #Absorption
                        lambAbs = 1/(ne_NaI*sigAbsNaI[r])
                        #NAbs = (round(DPshere2,2))*(1 - np.exp(-Ldetec*5.06e13/(lambAbs)))  # num of absorptions/year inside NaI crystals
                        NAbs = (DPshere2-NComps)*(1 - np.exp(-Ldetec*5.06e13/(lambAbs)))  # num of absorptions/year inside NaI crystals
                        if np.abs(NAbs) <=1e-12:
                            NAbs = 0
                            #print("sigAbs = %f, P = %f, E = %f" %(sigAbsNaI[r], 1 - np.exp(-Ldetec*5.06e13/(lambAbs)), E_z[r]))
                        # else:
                        #     print("NAbs = %.16f" %(NAbs))
                            #exit()
                        NAbs_tot[t][j][r] += NAbs

                        if np.sum(NDPs_tot[t][:j+1,r]) < np.sum(NAbs_tot[t][:j,r]):
                            print("j>0, AbsNaI", NAbs, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NAbs_tot[t][:j,r]), DPshere2, sigAbsNaI[r],t,j,E_z[r])
                            exit()
                            #NAbs_tot[t][j][r] -= np.abs(np.sum(NAbs_tot[t][:j,r])-np.sum(NDPs_tot[t][:j+1,r]))
                        
                        #Oscillation
                        Egamma=(E_z[r]+E_z[r-1])/2 #MeV
                        P=kappas**2*mDP**4*funcNaI(Egamma)*1e03*Ldetec*5.06e10/((mDP**2-(mgamma_NaI*1e03)**2)**2+Egamma**2*(funcNaI(Egamma)*1e03)**2)
                        if P<=1e-15:
                            P=0
                        if P>1:
                            P=1
                        NOsc = (DPshere2-NComps-NAbs)*P
                        if np.abs(NOsc) <=1e-12:
                            NOsc = 0
                        NOsc_tot[t][j][r] += NOsc 

                        if np.sum(NDPs_tot[t][:j+1,r]) < np.sum(NOsc_tot[t][:j,r]):
                            print("j>0, OscNaI", P, NOsc, DPshere2, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NOsc_tot[t][:j,r]),j, r)
                            exit()
                            #NOsc_tot[t][j][r] -= np.abs(np.sum(NOsc_tot[t][:j,r])-np.sum(NDPs_tot[t][:j+1,r]))

                        NDPs_res[t][r] -= (NComps+NAbs+NOsc)

                        #print(NComps, NAbs, NOsc)

                    else:

                        #Compton
                        lambComp = 1/(ne_NaI*sigCompNaI[r])
                        #NComps = round(NDPs_tot[t][j][r],1)*(1 - np.exp(-Ldetec*5.06e13/(lambComp))) # num of absorptions per year
                        NComps = (NDPs_res[t][r])*(1 - np.exp(-Ldetec*5.06e13/(lambComp))) # num of absorptions per year
                        if np.abs(NComps) <=1e-12:
                            NComps = 0
                            #print("sigComp = %f, P = %f, E = %f" %(sigCompNaI[r], 1 - np.exp(-Ldetec*5.06e13/(lambComp)), E_z[r]))
                        # else:
                        #     print("NComps = %.16f" %(NComps))
                            #exit()
                        
                        NComps_tot[t][j][r] += NComps

                        if NDPs_tot[t][j][r] < NComps_tot[t][j][r]:
                            print("j=0, CompNaI", NComps, lambComp, NDPs_tot[t][j][r], NComps_tot[t][j][r], np.exp(-L*5.06e13/(lambComp)), sigCompNaI[r],t,j,E_z[r])
                            exit()
                            NComps_tot[t][j][r] -= np.abs(NComps_tot[t][j][r]-NDPs_tot[t][j][r])    

                        #Absorption
                        lambAbs = 1/(ne_NaI*sigAbsNaI[r])
                        #NAbs = round(NDPs_tot[t][j][r],1)*(1 - np.exp(-Ldetec*5.06e13/(lambAbs)))#*1e-24 # num of absorptions per year
                        NAbs = (NDPs_res[t][r]-NComps)*(1 - np.exp(-Ldetec*5.06e13/(lambAbs)))#*1e-24 # num of absorptions per year
                        if np.abs(NAbs) <=1e-12:
                            NAbs = 0
                            #print("sigComp = %f, P = %f, E = %f" %(sigCompNaI[r], 1 - np.exp(-Ldetec*5.06e13/(lambComp)), E_z[r]))
                        #else:
                            #print("NAbs = %.16f" %(NAbs))
                            #exit()

                        NAbs_tot[t][j][r] += NAbs

                        if NDPs_tot[t][j][r] < NAbs_tot[t][j][r]:
                            print("j=0, AbsNaI", NAbs, lambAbs, NDPs_tot[t][j][r], NAbs_tot[t][j][r], np.exp(-L*5.06e13/(lambAbs)), sigAbsNaI[r],t,j,E_z[r])
                            exit()
                            NAbs_tot[t][j][r] -= np.abs(NAbs_tot[t][j][r]-NDPs_tot[t][j][r])

                        #Oscillation
                        NOsc = (NDPs_res[t][r]-NComps-NAbs)*P
                        if np.abs(NOsc) <=1e-12:
                            NOsc = 0

                        NOsc_tot[t][j][r] += NOsc 
                        
                        if NDPs_tot[t][j][r] < NOsc_tot[t][j][r]:
                            print("j=0, OscNaI", P, NOsc, DPshere2, j)
                            exit()
                            NOsc_tot[t][j][r] -= np.abs(NOsc_tot[t][j][r]-NDPs_tot[t][j][r])
                        NDPs_res[t][r] =- (NComps+NAbs+NOsc)

    return(NComps_tot,NAbs_tot,NOsc_tot)
    
def main():

    # Mass and epsilon values given by user. It makes this script to be parallelizable.
    mDP=float(sys.argv[1]) #MeV
    kappas=float(sys.argv[2])
    # mDP = 0.05963600 #MeV
    # kappas = 0.06493816316
    
    epsilon = 1e-4

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
    #mu=0.105658 #GeV  #muon mass
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
    
    # print(np.sum(Nmuabove100))
    # exit()

    #Building the mass attenuation function for each rock material
    (funcNa, en_Na)=InterpMassAtt("Na")
    (funcSi, en_Si)=InterpMassAtt("Si")
    (funcAl, en_Al)=InterpMassAtt("Al")
    (funcO, en_O)=InterpMassAtt("O")
    (funcFe, en_Fe)=InterpMassAtt("Fe")
    (funcCa, en_Ca)=InterpMassAtt("Ca")
    (funcK, en_K)=InterpMassAtt("K")
    (funcMg, en_Mg)=InterpMassAtt("Mg")
    (funcTi, en_Ti)=InterpMassAtt("Ti")
    (funcI, en_I)=InterpMassAtt("I")
    
    en1 = GatherEnergies(en_Na, en_Si)
    en2 = GatherEnergies(en1, en_Al)  
    en3 = GatherEnergies(en2, en_O)
    en4 = GatherEnergies(en3, en_Fe)
    en5 = GatherEnergies(en4, en_K)
    en6 = GatherEnergies(en5,en_Ca)
    en7 = GatherEnergies(en6,en_Ti)
    enfinal = GatherEnergies(en7, en_Mg)

    enfinal_NaI = GatherEnergies(en_Na, en_I)

    att_rock=np.zeros(len(enfinal))
    for i in range(len(enfinal)):
        att_rock[i]= 0.0283*funcNa(enfinal[i])+0.2772*funcSi(enfinal[i])+0.466*funcO(enfinal[i])+0.0813*funcAl(enfinal[i])
        +0.05*funcFe(enfinal[i])+0.0363*funcCa(enfinal[i])+0.0259*funcK(enfinal[i])+0.0209*funcMg(enfinal[i])+0.0148*funcTi(enfinal[i])
    funcrock=inter.interp1d(enfinal, att_rock*2.65*(1.973*10**(-14)), kind='linear', fill_value="extrapolate") # Mass attenuatioin for the rock in GeV, with energy given in MeV

    att_NaI=np.zeros(len(enfinal_NaI))  
    for i in range(len(enfinal_NaI)):
        att_NaI[i]= 0.5*funcNa(enfinal_NaI[i])+0.5*funcI(enfinal_NaI[i])
    funcNaI=inter.interp1d(enfinal_NaI, att_NaI*3.67*(1.973*10**(-14)), kind='linear', fill_value="extrapolate") # Mass attenuatioin for the rock in GeV, with energy given in MeV

    # Calculate how many DPs interact directly with the crystals.
    (NComps_tot, NAbs_tot, NOsc_tot) = CalcEvents(Ei, theta, E, mDP, x, kappas, Nmuabove100, Na, epsilon, funcrock, funcNaI)

    N = len(NComps_tot[0][0])
    NCompsNaI = np.zeros(N)
    NAbsNaI = np.zeros(N)
    NOscNaI = np.zeros(N)

    for i in range(0, N):
        NCompsNaI[i] = np.sum(NComps_tot[:,:,i])
        NAbsNaI[i] = np.sum(NAbs_tot[:,:,i])
        NOscNaI[i] = np.sum(NOsc_tot[:,:,i])

    if (mDP>1E-03):
        E_z = np.linspace(mDP, 0.2, N) #MeV
                    
    else:
        E_z = np.linspace(1E-03, 0.2, N) #MeV
        
    # with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/NDir100keV_int_m=%.8f_k=%.11f.txt' %(mDP,kappas), 'w') as archive:
    #     for i in range(0, N):
    #         archive.write('{:.11f}  {:.16f}  {:.16f}  {:.16f}\n'.format(E_z[i], NCompsNaI[i], NAbsNaI[i], NOscNaI[i]))
        
    with open('/home/dfreitas/DPnum_Comp/DirDetec/Spec100keV/NDir100keV_int_m=%.8f_k=%.11f.txt' %(mDP,kappas), 'w') as archive:
        for i in range(0, N):
            archive.write('{:.11f}  {:.16f}  {:.16f}  {:.16f}\n'.format(E_z[i], NCompsNaI[i], NAbsNaI[i], NOscNaI[i]))

main()