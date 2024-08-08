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
def diff_cross_section(eps, alpha, m_z, E_mu, m_mu, Z, x, A):
    '''
        Calculate the differential cross section for DP production by muon Bremsstrahlung
        using the IWW approximation
    '''
    m_z = m_z/1000
    chi_IWW = chi_calc(m_mu, E_mu, m_z, Z, A)
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
        Build the arrays of zenital angles and the muon ratios in this direction 
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
    #arq = open('/home/dfreitas/DPnum_Comp/Oscillation/Mass_att(Exmu)_%s.txt' %(element), 'r') 
    arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Oscillation/Mass_att(Exmu)_%s.txt' %(element), 'r') 
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

@numba.jit(nopython=True, nogil=True)
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

def CalcEvents(Ei, theta, E, mDP, x, kappas, Nmuabove100, Na, epsilon, funcmantle, funcore, funcNaI, R_earth):
    '''
    Function to compute the number NDPs_tot of muon Bremsstralung at the Earth
    producing DPs at each step j and angle t, for each coupling parameter and a
    DP mass mDP.
    '''
    alpha1=1/137
    me=0.511/1000 #electron mass (GeV)
    m_mu=105.7/1000 #muon mass (GeV)
    rho_mantle = 4 #mantle density (g/cm³) -> DOI: 10.1103/PhysRevLett.128.171801
    rho_core = 11  #core density (g/cm³) -> DOI: 10.1103/PhysRevLett.128.171801
    rho_NaI = 3.67 #NaI density (g/cm³)
    R_core = 3471000 #core's radius in (m)
    L_muon = np.diff(x)[0]/rho_mantle #step length inside the mantle for muon bremsstrahlung (cm)
    A_m = 24.5 #Mass number of the mantle (GeV) -> DOI: 10.1103/PhysRevLett.128.171801
    A_c = 54   #Mass number of the core (GeV) -> DOI: 10.1103/PhysRevLett.128.171801

    N=1000
    Z_mantle = 10.76
    Z_core = 26.2
    dsig_dx_m=np.zeros((len(Ei),N))
    sig_m=np.zeros((len(Ei),len(x)))
    ratio=np.zeros((len(Ei),N))

    #FIRST WE COMPUTE THE PROCESSES IN THE EARTH
    NDPs_tot = np.zeros((len(theta), N)) #Number of DPs generated in each angle and energy
    # NComps_res = np.zeros((len(theta), len(x), N)) #Number of Comptons in each angle, step and energy
    # NAbs_res = np.zeros((len(theta), len(x), N)) #Number of absorption processes in each angle, step and energy
    #NOsc_res = np.zeros((len(theta), 12000, N)) #Number of oscillations in each angle, step and energy in the Earth
    Ndec_res = np.zeros((len(theta), 12000, N)) #Number of decays in each angle, step and energy in Earth
    Ndec_tot = np.zeros((len(theta), 12000, N)) #Number of decays in each angle, step and energy in crystals

    q = 0.30282212088 #Elementary charge in natural units
    alpha = q**2/(4*np.pi) #Fine structure constant
    alphaline = (q*kappas)**2/(4*np.pi)

    ne_m=(rho_mantle/(A_m/Na))*(1.973*10**(-14))**3   #Electron number density for mantle (in GeV^3)
    ne_c=(rho_core/(A_c/Na))*(1.973*10**(-14))**3     #Electron number density for core (in GeV^3)
    ne_NaI=rho_NaI/(149.89/Na)*(1.973*10**(-14))**3   #Electron number density for NaI crystals (in GeV^3)
    mgamma_m=(4*np.pi*alpha1*ne_m/me)**0.5            #effective photon mass for mantle
    mgamma_c=(4*np.pi*alpha1*ne_c/me)**0.5            #effective photon mass for core

    #Compton-like scaterring Cross-Section
    # with open('/home/dfreitas/DPnum_Comp/Comp-xs2/CSxE_m=%f_k=%s00.txt' %(mDP, epsilon), 'r') as archive1:
    #         arqComp = archive1.read()

    with open('/home/davidc/Documents/Master\'s Analisys/Results Compton/Cross Section 2/CSxE_m=%f_k=%s00.txt' %(mDP, epsilon), 'r') as archive1:
            arqComp = archive1.read() 
        
    linesComp = arqComp.split("\n")
    linesComp = linesComp[:-1]

    for row in range(0, len(linesComp)):
        linesComp[row] = linesComp[row].split("   ")
        linesComp[row] = [float(x) for x in linesComp[row]]
    
    sigComparq_m = np.zeros(len(linesComp))
    sigComparq_c = np.zeros(len(linesComp))
    sigComparqNaI = np.zeros(len(linesComp))
    E_kComp = np.zeros(len(linesComp))

    for l in range(0, len(linesComp)):
        E_kComp[l] = linesComp[l][0]
        sigComparq_m[l] = (Z_mantle/64)*linesComp[l][1]*(1e-24)*(5.06e13)**2 #GeV⁻²/atom
        sigComparq_c[l] = (Z_core/64)*linesComp[l][1]*(1e-24)*(5.06e13)**2 #GeV⁻²/atom
        sigComparqNaI[l] = linesComp[l][1]*(1e-24)*(5.06e13)**2 #GeV⁻²/atom
    
    CScompfunc_m = inter.interp1d(E_kComp, sigComparq_m, kind = 'cubic', fill_value="extrapolate")
    CScompfunc_c = inter.interp1d(E_kComp, sigComparq_c, kind = 'cubic', fill_value="extrapolate")

    #Absorption Cross-Section
    with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/PeXSMantle.txt', 'r') as archivem:
            arqAbsm = archivem.read()
    with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/PeXSCore.txt', 'r') as archivec:
            arqAbsc = archivec.read() 
    # with open('/home/dfreitas/DPnum_Comp/Absorption/PeXSMantle.txt', 'r') as archivem:
    #         arqAbsm = archivem.read()
    # with open('/home/dfreitas/DPnum_Comp/Absorption/PeXSCore.txt', 'r') as archivec:
    #         arqAbsc = archivec.read()  

    linesAbsm = arqAbsm.split("\n")
    linesAbsc = arqAbsc.split("\n")
    linesAbsm = linesAbsm[:-1]
    linesAbsc = linesAbsc[:-1]
    
    for row in range(0, len(linesAbsm)):
        linesAbsm[row] = linesAbsm[row].split(" ")
        linesAbsm[row] = [float(x) for x in linesAbsm[row]]
    
    for row in range(0, len(linesAbsc)):
        linesAbsc[row] = linesAbsc[row].split(" ")
        linesAbsc[row] = [float(x) for x in linesAbsc[row]]

    sigAbs_m = np.zeros(len(linesAbsm))
    sigAbs_c = np.zeros(len(linesAbsc))
    E_kAbsm = np.zeros(len(linesAbsm))
    E_kAbsc = np.zeros(len(linesAbsc))

    for l in range(0, len(linesAbsm)):
        E_kAbsm[l] = linesAbsm[l][0]
        if mDP < E_kAbsm[l]:
            sigAbs_m[l] = (linesAbsm[l][1]*(5.06e13)**2/np.sqrt(1-mDP**2/E_kAbsm[l]**2))*(alphaline/alpha)*(A_m/Na) #GeV⁻²/atom (1/(Na/A)))
            
    for l in range(0, len(linesAbsc)):
            E_kAbsc[l] = linesAbsc[l][0]
            if mDP < E_kAbsc[l]:
                sigAbs_c[l] = (linesAbsc[l][1]*(5.06e13)**2/np.sqrt(1-mDP**2/E_kAbsc[l]**2))*(alphaline/alpha)*(A_c/Na) #GeV⁻²/atom (1/(Na/A)))

    CSabsfuncm = inter.interp1d(E_kAbsm, sigAbs_m, kind = 'linear', fill_value="extrapolate")
    CSabsfuncc = inter.interp1d(E_kAbsc, sigAbs_c, kind = 'linear', fill_value="extrapolate")
    CScompfuncNaI = inter.interp1d(E_kComp, sigComparqNaI, kind = 'cubic', fill_value="extrapolate")

    with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/PeXSNaI.txt', 'r') as archive3:
            arqAbs1 = archive3.read() 
    # with open('/home/dfreitas/DPnum_Comp/Absorption/PeXSNaI.txt', 'r') as archive3:
    #         arqAbs1 = archive3.read()

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

    #DP energies (GeV)
    if (mDP*1e-3>1E-06):
        E_z = np.linspace(mDP*1e-3, 0.1, N) #GeV
                    
    else:
        E_z = np.linspace(1E-06, 0.1, N) #GeV

    #FIRST, WE COMPUTE THE NUMBER OF DPs GENERATED BY MUON BREMSSTRAHLUNG IN EACH EARTH LAYER
    for i in range(len(Ei)-1):
        print("Bremsstrahlung", i)
        if (i>=900): continue
        for t in range(0,len(theta)):           
            for j in range(0,len(x)-1): 

                if ((14000-x[j]/rho_mantle/100)>=0.0):
                    if (E[i][j]<1*10**2): continue
                    for r in range(N):
                        ratio[i][r] = E_z[r]/(E[i][j]/1000)
            
                        if (r!=0):

                            #Bremsstrahlung
                            dsig_dx_m[i][r] = diff_cross_section(kappas, alpha1, mDP, E[i][j]/1000, m_mu, Z_mantle, (ratio[i][r]+ratio[i][r-1])/2, A_m)
                            integr=(dsig_dx_m[i][r]+dsig_dx_m[i][r-1])*(ratio[i][r]-ratio[i][r-1])/2
                            sig_m[i][j]=integr*(1.973*10**(-16))**2*10**4*10**24*10**12 #pb
                            xsec=sig_m[i][j]*10**(-12)*10**(-24) #cm^2
                            Nbremss=Nmuabove100[i][t][j]*rho_mantle*Na/A_m*L_muon*xsec*60*60*24*365.25 # num of DP's per year
                            NDPs_tot[t][r] += Nbremss

    totalDPs = np.sum(NDPs_tot)
            
    #NOW WE COMPUTE THE NUMBER OF STINGUISHED DPs INSIDE THE EARTH (DOESN'T WORK COMPUTING TOGETHER WITH THE BREMSSTRAHLUNG)
    sigCompm = (kappas/epsilon)**2*CScompfunc_m(1e03*E_z) #Compton cross-section
    sigAbsm = CSabsfuncm(1e03*E_z)  #Absorption cross section
    sigCompc = (kappas/epsilon)**2*CScompfunc_c(1e03*E_z) #Compton cross-section
    sigAbsc = CSabsfuncc(1e03*E_z)  #Absorption cross section
    sigAbsNaI = CSabsfuncNaI(E_z*1000)  #Absorption cross section in NaI
    sigCompNaI = (kappas/epsilon)**2*CScompfuncNaI(E_z) #Compton cross-section
    
    for t in range(0,len(theta)): #angles
        d = R_earth*np.sqrt(0.5*(1+np.cos(theta[t]))) #distance from the center to the sphere string
        print("events", t)
        if d >= R_core:
            xe = np.linspace(0,R_earth*np.sqrt(2*(1-np.cos(theta[t]))),12000) #Earth's path length (m)
            Le = np.diff(xe)[0]*100 #step length inside Earth (cm)
            print("\nmantle only")
            for j in range(0,len(xe)):  #steps inside rock
                if ((R_earth*np.sqrt(2*(1-np.cos(theta[t])))-xe[j])>=0.0):
                    
                    for r in range(0, N):

                        if r!= 0:

                            Egamma=(E_z[r]+E_z[r-1])/2 #GeV
                            lambAbsm = 1/(ne_m*sigAbsm[r])         #absorption length in mantle
                            lambAbsNaI = 1/(ne_NaI*sigAbsNaI[r])   #absorption length in NaI
                            lambCompm = 1/(ne_m*sigCompm[r])       #Compton length in mantle
                            lambCompNaI = 1/(ne_NaI*sigCompNaI[r]) #Compton length in NaI

                            if mDP/1000 > 2*me: #Computing the number of decays + events near the detector
                                gamma=kappas**2*alpha1*(mDP*1e-3)/3*(1+2*me**2/(mDP*1e-3)**2)*(1-4*me**2/(mDP*1e-3)**2)**0.5
                                Ldec = np.sqrt(E_z[r]**2/(mDP*1e-3)**2 - 1)/gamma  #decay length (GeV^{-1})
                                dist = (R_earth*np.sqrt(2*(1-np.cos(theta[t])))-xe[j])*100*5.06e13 #distance from step j to th detector ((GeV^{-1}))
                                Pdetec = np.exp(-dist*(1/Ldec+1/lambAbsm+1/lambCompm))*(1-np.exp(-4000*5.06e13*(1/Ldec+1/lambAbsNaI+1/lambCompNaI))) #detection probability (0.4 m detector) 
                                PEarth = 1 - np.exp(-Le*5.06e13*(1/lambAbsm+1/Ldec+1/lambCompm)) #Extinction probability in this step
                                Ndec_res[t][j][r] += PEarth*(round(NDPs_tot[t][r], 4)) #extinction number in this step
                                Ndec_tot[t][j][r] += Pdetec*(round(NDPs_tot[t][r], 4)) #detection number in this step via either direct
                                                                                        # interacion or decay  
                                NDPs_tot[t][r] -= (Ndec_res[t][j][r]+Ndec_tot[t][j][r])

                                if j>0:
                                    if round(NDPs_tot[t][r], 4) < round(np.sum(Ndec_res[t][:j,r]), 4):
                                        print(round(np.sum(Ndec_res[t][:j,r]), 4), round(NDPs_tot[t][r], 4))
                                        exit()
                                else:
                                    if round(NDPs_tot[t][r], 4) < round(Ndec_res[t][j][r], 4):
                                        print(round(Ndec_res[t][j][r], 4), round(NDPs_tot[t][r], 4))
                                        exit() 

                            else:
                                gamma = kappas**2*alpha1**4/(2**7*3**6*5**2*np.pi**3)*(mDP/1000)**9/me**8*(17/5+67/42*(mDP/1000)**2/me**2+128941/246960*(mDP/1000)**4/me**4)
                                Ldec = np.sqrt(E_z[r]**2/(mDP*1e-3)**2 - 1)/gamma
                                dist = (R_earth*np.sqrt(2*(1-np.cos(theta[t])))-xe[j])*100*5.06e13 #distance from step j to th detector ((GeV^{-1}))
                                P=kappas**2*(mDP*1e-03)**4*funcmantle(Egamma*1000)*Le*5.06e13/(((mDP*1e-03)**2-mgamma_m**2)**2+Egamma**2*funcmantle(Egamma*1000)**2)
                                Pdetec = np.exp(-dist*(1/Ldec+1/lambAbsm+1/lambCompm))*(1-np.exp(-4000*5.06e13*(1/Ldec+1/lambAbsNaI+1/lambCompNaI))) #detection probability (0.4 m detector) 
                                PEarth = 1 - np.exp(-Le*5.06e13*(1/lambAbsm+1/Ldec+1/lambCompm)) #Extinction probability in this step
                                NOsc = P*(round(NDPs_tot[t][r], 4))#Oscillation
                                Ndec_res[t][j][r] += PEarth*(round(NDPs_tot[t][r], 4)) + NOsc #extinction number in this step
                                Ndec_tot[t][j][r] += Pdetec*(round(NDPs_tot[t][r], 4))        #detection number in this step via either direct
                                    
                                NDPs_tot[t][r] -= (Ndec_res[t][j][r]+Ndec_tot[t][j][r])                                                     # interacion or decay 
                                if j>0:
                                    if round(NDPs_tot[t][r], 4) < round(np.sum(Ndec_res[t][:j,r]), 4):
                                        print(round(np.sum(Ndec_res[t][:j,r]), 4), round(NDPs_tot[t][r], 4))
                                        exit()
                                else:
                                    if round(NDPs_tot[t][r], 4) < round(Ndec_res[t][j][r], 4):
                                        print(round(Ndec_res[t][j][r], 4), round(NDPs_tot[t][r], 4))
                                        exit() 

                                # #Compton
                                # lambComp = 1/(ne_m*sigCompm[r])
                                # NComps = (round(DPshere, 4))*(1 - np.exp(-L_mantle*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                                # NComps_res[t][j][r] += NComps

                                # if round(np.sum(NDPs_tot[t][:j+1,r]), 4) < round(np.sum(NComps_res[t][:j,r]), 4):
                                #     # print("j>0, Comp", NComps, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NComps_res[t][:j,r]), np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), sigComp[r],t,j,E_z[r])
                                #     # exit()
                                #     NComps_res[t][j][r] -= NComps #REVER ESTA CONDIÇÃO!!!!!!!!!!!!!!!!

                                # #Absorption
                                # if mDP/1000 < 2*me:
                                #     lambAbs = 1/(ne_m*sigAbsm[r])
                                #     NAbs = (round(DPshere, 4))*(1 - np.exp(-L_mantle*5.06e13/(lambAbs)))  #num of Absorptions/year for this DP energy
                                #     NAbs_res[t][j][r] += NAbs
                                #     if round(np.sum(NDPs_tot[t][:j+1,r]), 4) < round(np.sum(NAbs_res[t][:j,r]), 4):
                                #         # print("j>0, Abs", NAbs, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NAbs_res[t][:j,r]), np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), sigAbs[r],t,j,E_z[r])
                                #         # exit()
                                #         NAbs_res[t][j][r] -= NAbs

                                # #Oscillation
                                # NOsc = round(DPshere,3)*P
                                # NOsc_res[t][j][r] += NOsc 

                                # if round(np.sum(NDPs_tot[t][:j+1,r]),3) < np.sum(NOsc_res[t][:j,r]):
                                #     # print("j>0, Osc", P, NOsc, round(DPshere,3), j)
                                #     # exit()
                                #     NOsc_res[t][j][r] -= NOsc
                            
                                # #Compton
                                # lambComp = 1/(ne_m*sigCompm[r])
                                # NComps = NDPs_tot[t][j][r]*(1 - np.exp(-L_mantle*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                                # NComps_res[t][j][r] += NComps

                                # if NComps_res[t][j][r] > NDPs_tot[t][j][r] or NComps > NDPs_tot[t][j][r]:
                                #     # print("j=0, Comp",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), NDPs_tot[t][j][r])
                                #     # exit()
                                #     NComps_res[t][j][r] -= NComps

                                # #Absorption
                                # if mDP/1000 < 2*me: 
                                #     lambAbs = 1/(ne_m*sigAbsm[r])
                                #     NAbs = NDPs_tot[t][j][r]*(1 - np.exp(-L_mantle*5.06e13/(lambAbs))) # num of absorptions per year for this DP energy
                                #     NAbs_res[t][j][r] += NAbs
                                #     if NAbs_res[t][j][r] > NDPs_tot[t][j][r] or NAbs > NDPs_tot[t][j][r]:
                                #         # print("j=0, Abs",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), NDPs_tot[t][j][r])
                                #         # exit()
                                #         NAbs_res[t][j][r] -= NAbs

                                # #Oscillation
                                # NOsc = NDPs_tot[t][j][r]*P
                                # NOsc_res[t][j][r] += NOsc

                                # if NDPs_tot[t][j][r] < NOsc:
                                #     # print("j=0, Osc", P, NOsc, DPshere)
                                #     # exit()
                                #     NOsc_res[t][j][r] -= NOsc
            
        else:
            L_c = 2*np.sqrt(R_core**2-R_earth**2/2*(1+np.cos(theta[t]))) #length of the core path (m)
            L_m = R_earth*np.sqrt(2*(1-np.cos(theta[t])))/2 - L_c/2 #length of each mantle path (m)
            x_m1 = np.linspace(0,L_m,4000)
            x_c = np.linspace(L_m, (L_m+L_c),4000)
            x_m2 = np.linspace((L_c+L_m),(L_c+2*L_m),4000)
            L_core = np.diff(x_c)[0]*100
            L_m1 = np.diff(x_m1)[0]*100
            L_m2 = np.diff(x_m2)[0]*100
            
            print("\n mantle 1")
            for j in range(0,len(x_m1)): #steps inside rock

                #Before recahing the nucleus
                if (L_m - x_m1[j] >= 0):
                    for r in range(0, N):
                        if r!= 0:
                            Egamma=(E_z[r]+E_z[r-1])/2 #GeV
                            lambAbsm = 1/(ne_m*sigAbsm[r])         #absorption length in mantle
                            lambAbsNaI = 1/(ne_NaI*sigAbsNaI[r])   #absorption length in NaI
                            lambCompm = 1/(ne_m*sigCompm[r])       #Compton length in mantle
                            lambCompNaI = 1/(ne_NaI*sigCompNaI[r]) #Compton length in NaI
                                
                            if mDP/1000 > 2*me: #Computing the number of decays + events near the detector
                                gamma=kappas**2*alpha1*(mDP*1e-3)/3*(1+2*me**2/(mDP*1e-3)**2)*(1-4*me**2/(mDP*1e-3)**2)**0.5
                                Ldec = np.sqrt(E_z[r]**2/(mDP*1e-3)**2 - 1)/gamma
                                dist = (R_earth*np.sqrt(2*(1-np.cos(theta[t])))-x_m1[j])*100*5.06e13 #distance from step j to th detector ((GeV^{-1}))
                                Pdetec = np.exp(-dist*(1/Ldec+1/lambAbsm+1/lambCompm))*(1-np.exp(-4000*5.06e13*(1/Ldec+1/lambAbsNaI+1/lambCompNaI))) #detection probability (0.4 m detector) 
                                PEarth = 1 - np.exp(-L_m1*5.06e13*(1/lambAbsm+1/Ldec+1/lambCompm)) #Extinction probability in this step
                                Ndec_res[t][j][r] += PEarth*(round(NDPs_tot[t][r], 4)) #extinction number in this step
                                Ndec_tot[t][j][r] += Pdetec*(round(NDPs_tot[t][r], 4)) #detection number in this step via either direct
                                                                                       # interacion or decay  
                                NDPs_tot[t][r] -= (Ndec_res[t][j][r]+Ndec_tot[t][j][r])
                                if j>0:
                                    if round(NDPs_tot[t][r], 4) < round(np.sum(Ndec_res[t][:j,r]), 4):
                                        print(round(np.sum(Ndec_res[t][:j,r]), 4), round(NDPs_tot[t][r], 4))
                                        exit()
                                else:
                                    if round(NDPs_tot[t][r], 4) < round(Ndec_res[t][j][r], 4):
                                        print(round(Ndec_res[t][j][r], 4), round(NDPs_tot[t][r], 4))
                                        exit() 
                            else:
                                gamma = kappas**2*alpha1**4/(2**7*3**6*5**2*np.pi**3)*(mDP/1000)**9/me**8*(17/5+67/42*(mDP/1000)**2/me**2+128941/246960*(mDP/1000)**4/me**4)
                                Ldec = np.sqrt(E_z[r]**2/(mDP*1e-3)**2 - 1)/gamma
                                dist = (R_earth*np.sqrt(2*(1-np.cos(theta[t])))-x_m1[j])*100*5.06e13 #distance from step j to th detector ((GeV^{-1}))
                                Pdetec = np.exp(-dist*(1/Ldec+1/lambAbsm+1/lambCompm))*(1-np.exp(-4000*5.06e13*(1/Ldec+1/lambAbsNaI+1/lambCompNaI))) #detection probability (0.4 m detector) 
                                PEarth = 1 - np.exp(-L_m1*5.06e13*(1/lambAbsm+1/Ldec+1/lambCompm)) #Extinction probability in this step
                                P=kappas**2*(mDP*1e-03)**4*funcmantle(Egamma*1000)*L_m1*5.06e13/(((mDP*1e-03)**2-mgamma_m**2)**2+Egamma**2*funcmantle(Egamma*1000)**2)
                                NOsc = P*(round(NDPs_tot[t][r], 4))#Oscillation
                                Ndec_res[t][j][r] += PEarth*(round(NDPs_tot[t][r], 4)) + NOsc #extinction number in this step
                                Ndec_tot[t][j][r] += Pdetec*(round(NDPs_tot[t][r], 4))        #detection number in this step via either direct
                                                                                              # interacion or decay 
                                NDPs_tot[t][r] -= (Ndec_res[t][j][r]+Ndec_tot[t][j][r])
                                if j>0:
                                    if round(NDPs_tot[t][r], 4) < round(np.sum(Ndec_res[t][:j,r]), 4):
                                        print(round(np.sum(Ndec_res[t][:j,r]), 4), round(NDPs_tot[t][r], 4))
                                        exit()
                                else:
                                    if round(NDPs_tot[t][r], 4) < round(Ndec_res[t][j][r], 4):
                                        print(round(Ndec_res[t][j][r], 4), round(NDPs_tot[t][r], 4))
                                        exit() 

                            # Egamma=(E_z[r]+E_z[r-1])/2 #GeV
                            # P=kappas**2*(mDP*1e-03)**4*funcmantle(Egamma*1000)*Le*5.06e13/(((mDP*1e-03)**2-mgamma_m**2)**2+Egamma**2*funcmantle(Egamma*1000)**2)

                            # if j>0:
                            #     DPshere = (np.sum(NDPs_tot[t][:j+1,r]) - np.sum(NComps_res[t][:j,r]) 
                            #                 - np.sum(NAbs_res[t][:j,r]) - np.sum(NOsc_res[t][:j,r]))
                            #     #Compton
                            #     lambComp = 1/(ne_m*sigCompm[r])
                            #     NComps = (round(DPshere, 4))*(1 - np.exp(-L_mantle*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                            #     NComps_res[t][j][r] += NComps

                            #     if round(np.sum(NDPs_tot[t][:j+1,r]), 4) < round(np.sum(NComps_res[t][:j,r]), 4):
                            #         # print("j>0, Comp", NComps, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NComps_res[t][:j,r]), np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), sigComp[r],t,j,E_z[r])
                            #         # exit()
                            #         NComps_res[t][j][r] -= NComps
                            
                            #     #Absorption
                            #     if mDP/1000 < 2*me: 
                            #         lambAbs = 1/(ne_m*sigAbsm[r])
                            #         NAbs = (round(DPshere, 4))*(1 - np.exp(-L_mantle*5.06e13/(lambAbs)))  #num of Absorptions/year for this DP energy
                            #         NAbs_res[t][j][r] += NAbs
                            #         if round(np.sum(NDPs_tot[t][:j+1,r]), 4) < round(np.sum(NAbs_res[t][:j,r]), 4):
                            #             # print("j>0, Abs", NAbs, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NAbs_res[t][:j,r]), np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), sigAbs[r],t,j,E_z[r])
                            #             # exit()
                            #             NAbs_res[t][j][r] -= NAbs

                            #     #Oscillation
                            #     NOsc = round(DPshere,3)*P
                            #     NOsc_res[t][j][r] += NOsc 

                            #     if round(np.sum(NDPs_tot[t][:j+1,r]),3) < np.sum(NOsc_res[t][:j,r]):
                            #         # print("j>0, Osc", P, NOsc, round(DPshere,3), j)
                            #         # exit()
                            #         NOsc_res[t][j][r] -= NOsc
                                        
                            # else:
                            #     #Compton
                            #     lambComp = 1/(ne_m*sigCompm[r])
                            #     NComps = NDPs_tot[t][j][r]*(1 - np.exp(-L_mantle*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                            #     NComps_res[t][j][r] += NComps

                            #     if NComps_res[t][j][r] > NDPs_tot[t][j][r] or NComps > NDPs_tot[t][j][r]:
                            #         # print("j=0, Comp",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), NDPs_tot[t][j][r])
                            #         # exit()
                            #         NComps_res[t][j][r] -= NComps

                            #     #Absorption
                            #     if mDP/1000 < 2*me: 
                            #         lambAbs = 1/(ne_m*sigAbsm[r])
                            #         NAbs = NDPs_tot[t][j][r]*(1 - np.exp(-L_mantle*5.06e13/(lambAbs))) # num of absorptions per year for this DP energy
                            #         NAbs_res[t][j][r] += NAbs
                            #         if NAbs_res[t][j][r] > NDPs_tot[t][j][r] or NAbs > NDPs_tot[t][j][r]:
                            #             # print("j=0, Abs",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), NDPs_tot[t][j][r])
                            #             # exit()
                            #             NAbs_res[t][j][r] -= NAbs

                            #     #Oscillation
                            #     NOsc = NDPs_tot[t][j][r]*P
                            #     NOsc_res[t][j][r] += NOsc

                            #     if NDPs_tot[t][j][r] < NOsc:
                            #         # print("j=0, Osc", P, NOsc, DPshere)
                            #         # exit()
                            #         NOsc_res[t][j][r] -= NOsc
            print("\n core")
            #In the nucleus
            for j in range(0,len(x_c)):
                if (L_m+L_c) - x_c[j] >= 0:
                    for r in range(0, N):
                        if r!= 0:
                        
                            Egamma=(E_z[r]+E_z[r-1])/2 #GeV
                            lambAbsc = 1/(ne_c*sigAbsc[r])         #absorption length in core
                            lambAbsNaI = 1/(ne_NaI*sigAbsNaI[r])   #absorption length in NaI
                            lambCompc = 1/(ne_c*sigCompc[r])       #Compton length in core
                            lambCompNaI = 1/(ne_NaI*sigCompNaI[r]) #Compton length in NaI
                                
                            if mDP/1000 > 2*me: #Computing the number of decays + events near the detector
                                gamma=kappas**2*alpha1*(mDP*1e-3)/3*(1+2*me**2/(mDP*1e-3)**2)*(1-4*me**2/(mDP*1e-3)**2)**0.5
                                Ldec = np.sqrt(E_z[r]**2/(mDP*1e-3)**2 - 1)/gamma
                                dist = (R_earth*np.sqrt(2*(1-np.cos(theta[t])))-x_c[j])*100*5.06e13 #distance from step j to th detector ((GeV^{-1}))
                                Pdetec = np.exp(-dist*(1/Ldec+1/lambAbsc+1/lambCompc))*(1-np.exp(-4000*5.06e13*(1/Ldec+1/lambAbsNaI+1/lambCompNaI))) #detection probability (0.4 m detector) 
                                PEarth = 1 - np.exp(-L_core*5.06e13*(1/lambAbsm+1/Ldec+1/lambCompc)) #Extinction probability in this step
                                Ndec_res[t][j+len(x_m1)][r] += PEarth*(round(NDPs_tot[t][r], 4))  #extinction number in this step
                                Ndec_tot[t][j+len(x_m1)][r] += Pdetec*(round(NDPs_tot[t][r], 4))  #detection number in this step via either direct
                                                                                                  # interacion or decay  

                                NDPs_tot[t][r] -= (Ndec_res[t][j][r]+Ndec_tot[t][j][r])
                                if j>0:
                                    if round(NDPs_tot[t][r], 4) < round(np.sum(Ndec_res[t][:j,r]), 4):
                                        print(round(np.sum(Ndec_res[t][:j,r]), 4), round(NDPs_tot[t][r], 4))
                                        exit()
                                else:
                                    if round(NDPs_tot[t][r], 4) < round(Ndec_res[t][j][r], 4):
                                        print(round(Ndec_res[t][j][r], 4), round(NDPs_tot[t][r], 4))
                                        exit() 

                            else:
                                gamma = kappas**2*alpha1**4/(2**7*3**6*5**2*np.pi**3)*(mDP/1000)**9/me**8*(17/5+67/42*(mDP/1000)**2/me**2+128941/246960*(mDP/1000)**4/me**4)
                                Ldec = np.sqrt(E_z[r]**2/(mDP*1e-3)**2 - 1)/gamma
                                dist = (R_earth*np.sqrt(2*(1-np.cos(theta[t])))-x_c[j])*100*5.06e13 #distance from step j to th detector ((GeV^{-1}))
                                Pdetec = np.exp(-dist*(1/Ldec+1/lambAbsc+1/lambCompc))*(1-np.exp(-4000*5.06e13*(1/Ldec+1/lambAbsNaI+1/lambCompNaI))) #detection probability (0.4 m detector) 
                                PEarth = 1 - np.exp(-L_core*5.06e13*(1/lambAbsc+1/Ldec+1/lambCompc)) #Extinction probability in this step
                                P=kappas**2*(mDP*1e-03)**4*funcore(Egamma*1000)*L_core*5.06e13/(((mDP*1e-03)**2-mgamma_c**2)**2+Egamma**2*funcore(Egamma*1000)**2)
                                NOsc = P*(round(NDPs_tot[t][r], 4)) #Oscillation
                                Ndec_res[t][j+len(x_m1)][r] += PEarth*(round(NDPs_tot[t][r], 4)) + NOsc #extinction number in this step
                                Ndec_tot[t][j+len(x_m1)][r] += Pdetec*(round(NDPs_tot[t][r], 4))        #detection number in this step via either direct
                                                                                                        # interacion or decay 
                            
                                NDPs_tot[t][r] -= (Ndec_res[t][j][r]+Ndec_tot[t][j][r])
                                if j>0:
                                    if round(NDPs_tot[t][r], 4) < round(np.sum(Ndec_res[t][:j,r]), 4):
                                        print(round(np.sum(Ndec_res[t][:j,r]), 4), round(NDPs_tot[t][r], 4))
                                        exit()
                                else:
                                    if round(NDPs_tot[t][r], 4) < round(Ndec_res[t][j][r], 4):
                                        print(round(Ndec_res[t][j][r], 4), round(NDPs_tot[t][r], 4))
                                        exit() 

                                # DPshere = (np.sum(NDPs_tot[t][:j+1,r]) - np.sum(NComps_res[t][:j,r]) 
                                #             - np.sum(NAbs_res[t][:j,r]) - np.sum(NOsc_res[t][:j,r]))
                                # #Compton
                                # lambComp = 1/(ne_c*sigCompc[r])
                                # NComps = (round(DPshere, 4))*(1 - np.exp(-L_core*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                                # NComps_res[t][j][r] += NComps

                                # if round(np.sum(NDPs_tot[t][:j+1,r]), 4) < round(np.sum(NComps_res[t][:j,r]), 4):
                                #     # print("j>0, Comp", NComps, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NComps_res[t][:j,r]), np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), sigComp[r],t,j,E_z[r])
                                #     # exit()
                                #     NComps_res[t][j][r] -= NComps
                            
                                # #Absorption
                                # if mDP/1000 < 2*me: 
                                #     lambAbs = 1/(ne_c*sigAbsc[r])
                                #     NAbs = (round(DPshere, 4))*(1 - np.exp(-L_core*5.06e13/(lambAbs)))  #num of Absorptions/year for this DP energy
                                #     NAbs_res[t][j][r] += NAbs
                                #     if round(np.sum(NDPs_tot[t][:j+1,r]), 4) < round(np.sum(NAbs_res[t][:j,r]), 4):
                                #         # print("j>0, Abs", NAbs, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NAbs_res[t][:j,r]), np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), sigAbs[r],t,j,E_z[r])
                                #         # exit()
                                #         NAbs_res[t][j][r] -= NAbs

                                # #Oscillation
                                # NOsc = round(DPshere,3)*P
                                # NOsc_res[t][j][r] += NOsc 

                                # if round(np.sum(NDPs_tot[t][:j+1,r]),3) < np.sum(NOsc_res[t][:j,r]):
                                #     # print("j>0, Osc", P, NOsc, round(DPshere,3), j)
                                #     # exit()
                                #     NOsc_res[t][j][r] -= NOsc
                            
                                # #Compton
                                # lambComp = 1/(ne_c*sigCompc[r])
                                # NComps = NDPs_tot[t][j][r]*(1 - np.exp(-L_core*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                                # NComps_res[t][j][r] += NComps

                                # if NComps_res[t][j][r] > NDPs_tot[t][j][r] or NComps > NDPs_tot[t][j][r]:
                                #     # print("j=0, Comp",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), NDPs_tot[t][j][r])
                                #     # exit()
                                #     NComps_res[t][j][r] -= NComps

                                # #Absorption
                                # if mDP/1000 < 2*me:
                                #     lambAbs = 1/(ne_c*sigAbsc[r])
                                #     NAbs = NDPs_tot[t][j][r]*(1 - np.exp(-L_core*5.06e13/(lambAbs))) # num of absorptions per year for this DP energy
                                #     NAbs_res[t][j][r] += NAbs
                                #     if NAbs_res[t][j][r] > NDPs_tot[t][j][r] or NAbs > NDPs_tot[t][j][r]:
                                #         # print("j=0, Abs",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), NDPs_tot[t][j][r])
                                #         # exit()
                                #         NAbs_res[t][j][r] -= NAbs

                                # #Oscillation
                                # NOsc = NDPs_tot[t][j][r]*P
                                # NOsc_res[t][j][r] += NOsc

                                # if NDPs_tot[t][j][r] < NOsc:
                                #     # print("j=0, Osc", P, NOsc, DPshere)
                                #     # exit()
                                #     NOsc_res[t][j][r] -= NOsc

            print("\n mantle 2")
            #After coming out of the nucleus
            for j in range(0, len(x_m2)):
                if R_earth*np.sqrt(2*(1-np.cos(theta[t])))-x_m2[j]>=0:
                    for r in range(0, N):
                        if r!= 0:

                            Egamma=(E_z[r]+E_z[r-1])/2 #GeV
                            lambAbsm = 1/(ne_m*sigAbsm[r])         #absorption length in mantle
                            lambAbsNaI = 1/(ne_NaI*sigAbsNaI[r])   #absorption length in NaI
                            lambCompm = 1/(ne_m*sigCompm[r])       #Compton length in mantle
                            lambCompNaI = 1/(ne_NaI*sigCompNaI[r]) #Compton length in NaI

                            if mDP/1000 > 2*me: #Computing the number of decays + events near the detector
                                gamma=kappas**2*alpha1*(mDP*1e-3)/3*(1+2*me**2/(mDP*1e-3)**2)*(1-4*me**2/(mDP*1e-3)**2)**0.5
                                Ldec = np.sqrt(E_z[r]**2/(mDP*1e-3)**2 - 1)/gamma
                                dist = (R_earth*np.sqrt(2*(1-np.cos(theta[t])))-x_m2[j])*100*5.06e13 #distance from step j to th detector ((GeV^{-1}))
                                Pdetec = np.exp(-dist*(1/Ldec+1/lambAbsm+1/lambCompm))*(1-np.exp(-4000*5.06e13*(1/Ldec+1/lambAbsNaI+1/lambCompNaI))) #detection probability (0.4 m detector) 
                                PEarth = 1 - np.exp(-L_m2*5.06e13*(1/lambAbsm+1/Ldec+1/lambCompm)) #Extinction probability in this step
                                Ndec_res[t][j+len(x_m1)+len(x_c)][r] += PEarth*(round(NDPs_tot[t][r], 4)) #extinction number in this step
                                Ndec_tot[t][j+len(x_m1)+len(x_c)][r] += Pdetec*(round(NDPs_tot[t][r], 4)) #detection number in this step via either direct
                                                                                                          # interacion or decay  

                                NDPs_tot[t][r] -= (Ndec_res[t][j][r]+Ndec_tot[t][j][r])
                                if j>0:
                                    if round(NDPs_tot[t][r], 4) < round(np.sum(Ndec_res[t][:j,r]), 4):
                                        print(round(np.sum(Ndec_res[t][:j,r]), 4), round(NDPs_tot[t][r], 4))
                                        exit()
                                else:
                                    if round(NDPs_tot[t][r], 4) < round(Ndec_res[t][j][r], 4):
                                        print(round(Ndec_res[t][j][r], 4), round(NDPs_tot[t][r], 4))
                                        exit() 
                            else:
                                gamma = kappas**2*alpha1**4/(2**7*3**6*5**2*np.pi**3)*(mDP/1000)**9/me**8*(17/5+67/42*(mDP/1000)**2/me**2+128941/246960*(mDP/1000)**4/me**4)
                                Ldec = np.sqrt(E_z[r]**2/(mDP*1e-3)**2 - 1)/gamma
                                dist = (R_earth*np.sqrt(2*(1-np.cos(theta[t])))-x_m2[j])*100*5.06e13 #distance from step j to th detector ((GeV^{-1}))
                                Pdetec = np.exp(-dist*(1/Ldec+1/lambAbsm+1/lambCompm))*(1-np.exp(-4000*5.06e13*(1/Ldec+1/lambAbsNaI+1/lambCompNaI))) #detection probability (0.4 m detector) 
                                PEarth = 1 - np.exp(-L_m2*5.06e13*(1/lambAbsm+1/Ldec+1/lambCompm)) #Extinction probability in this step
                                P=kappas**2*(mDP*1e-03)**4*funcmantle(Egamma*1000)*L_m2*5.06e13/(((mDP*1e-03)**2-mgamma_m**2)**2+Egamma**2*funcmantle(Egamma*1000)**2)
                                NOsc = P*(round(NDPs_tot[t][r], 4))#Oscillation
                                Ndec_res[t][j+len(x_m1)+len(x_c)][r] += PEarth*(round(NDPs_tot[t][r], 4)) + NOsc #extinction number in this step
                                Ndec_tot[t][j+len(x_m1)+len(x_c)][r] += Pdetec*(round(NDPs_tot[t][r], 4))        #detection number in this step via either direct
                                                                                                                 # interacion or decay 
                                    
                                NDPs_tot[t][r] -= (Ndec_res[t][j][r]+Ndec_tot[t][j][r])
                                if j>0:
                                    if round(NDPs_tot[t][r], 4) < round(np.sum(Ndec_res[t][:j,r]), 4):
                                        print(round(np.sum(Ndec_res[t][:j,r]), 4), round(NDPs_tot[t][r], 4))
                                        exit()
                                else:
                                    if round(NDPs_tot[t][r], 4) < round(Ndec_res[t][j][r], 4):
                                        print(round(Ndec_res[t][j][r], 4), round(NDPs_tot[t][r], 4))
                                        exit() 

                                # DPshere = (np.sum(NDPs_tot[t][:j+1,r]) - np.sum(NComps_res[t][:j,r]) 
                                #             - np.sum(NAbs_res[t][:j,r]) - np.sum(NOsc_res[t][:j,r]))
                                # #Compton
                                # lambComp = 1/(ne_m*sigCompm[r])
                                # NComps = (round(DPshere, 4))*(1 - np.exp(-L_mantle*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                                # NComps_res[t][j][r] += NComps

                                # if round(np.sum(NDPs_tot[t][:j+1,r]), 4) < round(np.sum(NComps_res[t][:j,r]), 4):
                                #     # print("j>0, Comp", NComps, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NComps_res[t][:j,r]), np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), sigComp[r],t,j,E_z[r])
                                #     # exit()
                                #     NComps_res[t][j][r] -= NComps
                            
                                # #Absorption
                                # if mDP/1000 < 2*me: 
                                #     lambAbs = 1/(ne_m*sigAbsm[r])
                                #     NAbs = (round(DPshere, 4))*(1 - np.exp(-L_mantle*5.06e13/(lambAbs)))  #num of Absorptions/year for this DP energy
                                #     NAbs_res[t][j][r] += NAbs
                                #     if round(np.sum(NDPs_tot[t][:j+1,r]), 4) < round(np.sum(NAbs_res[t][:j,r]), 4):
                                #         # print("j>0, Abs", NAbs, lambAbs, np.sum(NDPs_tot[t][:j+1,r]), np.sum(NAbs_res[t][:j,r]), np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), sigAbs[r],t,j,E_z[r])
                                #         # exit()
                                #         NAbs_res[t][j][r] -= NAbs

                                # #Oscillation
                                # NOsc = round(DPshere,3)*P
                                # NOsc_res[t][j][r] += NOsc 

                                # if round(np.sum(NDPs_tot[t][:j+1,r]),3) < np.sum(NOsc_res[t][:j,r]):
                                #     # print("j>0, Osc", P, NOsc, round(DPshere,3), j)
                                #     # exit()
                                #     NOsc_res[t][j][r] -= NOsc
                            
                                # #Compton
                                # lambComp = 1/(ne_m*sigCompm[r])
                                # NComps = NDPs_tot[t][j][r]*(1 - np.exp(-L_mantle*5.06e13/(lambComp)))  #num of Comptons/year for this DP energy
                                # NComps_res[t][j][r] += NComps

                                # if NComps_res[t][j][r] > NDPs_tot[t][j][r] or NComps > NDPs_tot[t][j][r]:
                                #     # print("j=0, Comp",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), NDPs_tot[t][j][r])
                                #     # exit()
                                #     NComps_res[t][j][r] -= NComps

                                # #Absorption
                                # if mDP/1000 < 2*me: 
                                #     lambAbs = 1/(ne_m*sigAbsm[r])
                                #     NAbs = NDPs_tot[t][j][r]*(1 - np.exp(-L_mantle*5.06e13/(lambAbs))) # num of absorptions per year for this DP energy
                                #     NAbs_res[t][j][r] += NAbs
                                #     if NAbs_res[t][j][r] > NDPs_tot[t][j][r] or NAbs > NDPs_tot[t][j][r]:
                                #         # print("j=0, Abs",NDPs_tot[t][j][r]*np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), NAbs_res[t][j][r], np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), NDPs_tot[t][j][r])
                                #         # exit()
                                #         NAbs_res[t][j][r] -= NAbs

                                # #Oscillation
                                # NOsc = NDPs_tot[t][j][r]*P
                                # NOsc_res[t][j][r] += NOsc

                                # if NDPs_tot[t][j][r] < NOsc:
                                #     # print("j=0, Osc", P, NOsc, DPshere)
                                #     # exit()
                                #     NOsc_res[t][j][r] -= NOsc

########################################################################################################

    for t in range(0,len(theta)): #angles
        for r in range(0, N):
            if r!= 0:
                NDPs_tot[t][r] = round(NDPs_tot[t][r]-np.sum(Ndec_res[t][:,r]), 1)

                        # if NDPs_tot[t][j][r] < 0:
                        #     print(NDPs_tot[t][j][r])
                        #     NDPs_tot[t][j][r] = 0
    #NDPs_tot = round(np.sum(NDPs_tot), 4)

    #NOW WE COMPUTE THE INTERACTIONS WITH THE CRYSTAL
    if mDP/1000 < 2*me:
        detectheta = np.arcsin(np.sin(theta)/(np.sqrt(2*(1-np.cos(theta)))))
        xdetect = np.linspace(0,0.4/np.cos(detectheta[len(detectheta)-1]),500) #Steps inside the NaI crystals (m)
        # NComps_tot = np.zeros((len(theta), len(xdetect), N)) #Number of Comptons in each angle, step and energy
        # NAbs_tot = np.zeros((len(theta), len(xdetect), N)) #Number of absorptions in each angle, step and energy
        NOsc_tot = np.zeros((len(detectheta), len(xdetect), N)) #Number of oscillations in each angle, step and energy
        Ldetec = np.diff(xdetect)[0] 
        mgamma_NaI=(4*np.pi*alpha1*ne_NaI/me)**0.5    #effective photon mass (GeV)

        if (mDP>1E-03):
            E_z = np.linspace(mDP, 3, N) #MeV
                        
        else:
            E_z = np.linspace(1E-03, 3, N) #MeV

        for t in range(0,len(theta)):
            print("detection", t)
            for j in range(0,len(xdetect)-1):
                if ((0.4/np.cos(theta[t])-xdetect[j])>=0.0):

                    for r in range(N):

                        #if j>0:
                            #DPshere2 = (NDPs_tot[t][r] - np.sum(NOsc_tot[t][:j,r]))
        #                     #Compton
        #                     lambComp = 1/(ne_NaI*sigCompNaI[r])
        #                     NComps = (round(DPshere2,2))*(1 - np.exp(-Ldetec*5.06e13/(lambComp))) # num of Compton/year inside NaI crystals
        #                     NComps_tot[t][j][r] += NComps

        #                     if round(np.sum(NDPs_tot[t][:j+1,r]), 2) < round(np.sum(NComps_tot[t][:j,r]), 2):
        #                         # print("j>0, CompNaI", NComps, lambAbs, round(np.sum(NDPs_tot[t][:j+1,r]), 2), round(np.sum(NComps_tot[t][:j,r]),2), round(DPshere2,2), sigCompNaI[r],t,j,E_z[r])
        #                         # exit()
        #                         NComps_tot[t][j][r] -= NComps

        #                     #Absorption
        #                     lambAbs = 1/(ne_NaI*sigAbsNaI[r])
        #                     NAbs = (round(DPshere2,2))*(1 - np.exp(-Ldetec*5.06e13/(lambAbs)))  # num of absorptions/year inside NaI crystals
        #                     NAbs_tot[t][j][r] += NAbs

        #                     if round(np.sum(NDPs_tot[t][:j+1,r]), 2) < round(np.sum(NAbs_tot[t][:j,r]), 2):
        #                         # print("j>0, AbsNaI", NAbs, lambAbs, round(np.sum(NDPs_tot[t][:j+1,r]),2), round(np.sum(NAbs_tot[t][:j,r]),2), round(DPshere2,3), sigAbsNaI[r],t,j,E_z[r])
        #                         # exit()
        #                         NAbs_tot[t][j][r] -= NAbs

                        #Oscillation
                        Egamma=(E_z[r]+E_z[r-1])/2 #MeV
                        P=kappas**2*mDP**4*funcNaI(Egamma)*1e03*Ldetec*5.06e10/((mDP**2-(mgamma_NaI*1e03)**2)**2+Egamma**2*(funcNaI(Egamma)*1e03)**2)
                        NOsc = round(NDPs_tot[t][r],3)*P
                        NOsc_tot[t][j][r] += NOsc 
                        NDPs_tot[t][r] -= NOsc

                        if j>0:
                            if round(NDPs_tot[t][r],3) < np.sum(NOsc_tot[t][:j,r]):
                                # print("j>0, Osc", P, NOsc, round(DPshere,3), j)
                                # exit()
                                NOsc_tot[t][j][r] -= NOsc
                        else:
                            if round(NDPs_tot[t][r],3) < NOsc_tot[t][j][r]:
                                # print("j>0, Osc", P, NOsc, round(DPshere,3), j)
                                # exit()
                                NOsc_tot[t][j][r] -= NOsc
                        #else:
        #                     #Compton
        #                     lambComp = 1/(ne_NaI*sigCompNaI[r])
        #                     NComps = round(NDPs_tot[t][j][r],1)*(1 - np.exp(-Ldetec*5.06e13/(lambComp))) # num of absorptions per year
        #                     NComps_tot[t][j][r] += NComps

        #                     if round(NDPs_tot[t][j][r], 2) < round(NComps_tot[t][j][r], 2):
        #                     #     print("j=0, CompNaI", NComps, lambComp, round(NDPs_tot[t][j][r],2), round(NComps_tot[t][j][r],2), np.exp(-L*5.06e13/(lambComp*np.cos(theta[t]))), sigCompNaI[r],t,j,E_z[r])
        #                     #     exit()
        #                         NComps_tot[t][j][r] -= NComps    

        #                     #Absorption
        #                     lambAbs = 1/(ne_NaI*sigAbsNaI[r])
        #                     NAbs = round(NDPs_tot[t][j][r],1)*(1 - np.exp(-Ldetec*5.06e13/(lambAbs)))#*1e-24 # num of absorptions per year
        #                     NAbs_tot[t][j][r] += NAbs

        #                     if round(NDPs_tot[t][j][r], 1) < round(NAbs_tot[t][j][r], 1):
        #                     #     print("j=0, AbsNaI", NAbs, lambAbs, round(NDPs_tot[t][j][r],1), round(NAbs_tot[t][j][r],1), np.exp(-L*5.06e13/(lambAbs*np.cos(theta[t]))), sigAbsNaI[r],t,j,E_z[r])
        #                     #     exit()
        #                         NAbs_tot[t][j][r] -= NAbs
                            
                            # #Oscillation
                            # Egamma=(E_z[r]+E_z[r-1])/2 #MeV
                            # P=kappas**2*mDP**4*funcNaI(Egamma)*Ldetec*5.06e10/((mDP**2-mgamma_NaI**2)**2+Egamma**2*funcNaI(Egamma)**2)
                            # NOsc = round(NDPs_tot[t][r],1)*P
                            # NOsc_tot[t][j][r] += NOsc 

                            # if round(np.sum(NDPs_tot[t][:,r]),3) < np.sum(NOsc_tot[t][:j,r]):
                            #     # print("j>0, Osc", P, NOsc, round(DPshere,3), j)
                            #     # exit()
                            #     NOsc_tot[t][j][r] -= NOsc

    return(NDPs_tot, Ndec_tot, Ndec_res, totalDPs, NOsc_tot)
    
def main():

    # Mass and epsilon values given by user. It makes this script to be parallelizable.
    # mDP=float(sys.argv[1]) #MeV
    # kappas=float(sys.argv[2])
    
    mDP = 0.000121 #MeV
    kappas = 0.00133352143
    epsilon = 1e-4
    R_earth = 6371000   #Earth's radius in meters

    # Import vertical muon flux from around 0.8 GeV to few TeV
    arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonflux.tsv', 'r')
    #arq = open('/home/dfreitas/DPnum_Comp/Muons/muonflux.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    p=np.zeros(N1)    #muon momentum
    y=np.zeros(N1)
    #E=np.zeros(N1)    #muon energy
    mu=0.105658 #GeV  #muon mass
    Na=6.02*10**23    #Avogadro number

    for i in range(N1):
        p[i]=texto1[i][0]
        y[i]=texto1[i][1]
        #E[i]=(p[i]**2+mu**2)**0.5 # Since p>>mu, E~p
   
    # Interpolate muon flux in units of muons/cm^2/s/sr
    func=inter.interp1d(p,y/p**3)

    # Import muon stopping power in rock
    arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muon_stoppower.tsv', 'r')
    #arq = open('/home/dfreitas/DPnum_Comp/Muons/muon_stoppower.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    dEdx=np.zeros(N1) #Energy variation rate for the muon in each step
    Em=np.zeros(N1)   

    Build_general(texto1, Em, dEdx, N1)

    funcdEdx=inter.interp1d(Em,dEdx)

    # Import muon flux variation with the zenith angle for 1 GeV muons
    arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonangle_1GeV.tsv', 'r')
    #arq = open('/home/dfreitas/DPnum_Comp/Muons/muonangle_1GeV.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    zenite_1GeV=np.zeros(N1)
    ratio_1GeV=np.zeros(N1)

    Build_zenital(texto1, zenite_1GeV, ratio_1GeV, N1)

    func1GeV=inter.interp1d(zenite_1GeV,ratio_1GeV)

    # Import muon flux variation with the zenith angle for 10 GeV muons
    arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonangle_10GeV.tsv', 'r')
    #arq = open('/home/dfreitas/DPnum_Comp/Muons/muonangle_10GeV.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    zenite_10GeV=np.zeros(N1)
    ratio_10GeV=np.zeros(N1)

    Build_zenital(texto1, zenite_10GeV, ratio_10GeV, N1)

    func10GeV=inter.interp1d(zenite_10GeV,ratio_10GeV)

    # Import muon flux variation with the zenith angle for 100 GeV muons
    arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonangle_100GeV.tsv', 'r')
    #arq = open('/home/dfreitas/DPnum_Comp/Muons/muonangle_100GeV.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    zenite_100GeV=np.zeros(N1)
    ratio_100GeV=np.zeros(N1)

    Build_zenital(texto1, zenite_100GeV, ratio_100GeV, N1)

    func100GeV=inter.interp1d(zenite_100GeV,ratio_100GeV)

    # Import muon flux variation with the zenith angle for 1 TeV muons
    arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/muonangle_1000GeV.tsv', 'r')
    #arq = open('/home/dfreitas/DPnum_Comp/Muons/muonangle_1000GeV.tsv', 'r')
    texto1 = arq.readlines()
    arq.close()

    texto1 = ExtractData(texto1)

    N1=len(texto1)
    zenite_1000GeV=np.zeros(N1)
    ratio_1000GeV=np.zeros(N1)

    Build_zenital(texto1, zenite_1000GeV, ratio_1000GeV, N1)

    func1000GeV=inter.interp1d(zenite_1000GeV,ratio_1000GeV)

    # Interpolate muon flux variation with zenith for different muon energies
    theta_earth = np.linspace(3.0*np.pi/180,177*np.pi/180,60) # Earth's theta angle in spherical coordinates, considering 3 degrees resolution
    
    #print(theta*180/np.pi)
    en=np.logspace(0,3,4)
    interen=np.zeros(4)
    Ei=np.logspace(3,6.69897,1000) #Muon energies from 1 GeV to few TeV (in MeV)
    #print((Ei[0]-np.diff(Ei)[0]/2)/1000)
    ratio=np.zeros((len(theta_earth),len(Ei)))

    for t in range(len(theta_earth)):
        dvec = np.array([-R_earth*np.sin(theta_earth[t]),0,R_earth*(1-np.cos(theta_earth[t]))])
        rhat = np.array([np.sin(theta_earth[t]),0,np.cos(theta_earth[t])])
        zenite = np.arccos(np.dot(-dvec,rhat)/(R_earth*np.sqrt(2*(1-np.cos(theta_earth[t])))))
        interen[0]=func1GeV(zenite)
        interen[1]=func10GeV(zenite)
        interen[2]=func100GeV(zenite)
        interen[3]=func1000GeV(zenite)
        funcen=inter.interp1d(en,interen)
        for i in range(len(Ei)):
            if (Ei[i]<=1E+6):
                ratio[t][i]=funcen(Ei[i]/1000)
            else:
                ratio[t][i]=2.0*funcen(1E+6/1000) # Since the muon flux variation with zenith is only given to energies up to 1 TeV, consider that for energies above
                                                  # 1 TeV the muon flux is two times the muon flux for 1 TeV.

    #x=np.linspace(0,(700/np.cos(theta[14]))/700*1855*100,5000) # Define step length in cm used in the muon attenuation and DP generated by Bremsstrahlung in the rock.
    x=np.linspace(0,14000*2.65*100,6000) # Define step length in cm used in the muon attenuation and DP generated by Bremsstrahlung in the rock.

    # Calculate muon attenuation when travelling in the rock
    E=np.zeros((len(Ei),len(x)))
    for i in range(len(Ei)):
        
        for j in range(len(x)-1):
            if x[j]/2.65/100 <= 14000: #Rock density = 2.65 g/cm³ and muons penetrate till' 14000 m max
                if (j==0):
                    E[i][j]=Ei[i]

                if (E[i][j]<=1.2): continue
                E[i][j+1]=E[i][j]-funcdEdx(E[i][j])*(x[j+1]-x[j]) 

    Acosine=(2*100)**2 #cm^2 COSINE-100 LS Area
    Nmu=np.zeros((len(theta_earth),len(Ei)-1))

    # Check how many muons had a trajectory that would cross COSINE LS detector
    for i in range(len(Ei)-1):

        for t in range(len(theta_earth)): 
            I = quad(func, (Ei[i]-np.diff(Ei)[i]/2)/1000, (Ei[i]+np.diff(Ei)[i]/2)/1000)
            sr = Acosine/(2*(R_earth*100)**2*(1-np.cos(theta_earth[t]))) #Omega = Area/distance²
            Nmu[t][i] = ratio[t][i]*I[0]*Acosine*sr #Muons/s

    # Check how many muons there are in each zenith angle, depth in rock, and with energy E[i][j]. 
    # It is required that the muon energy has to be greater than 100 MeV.
    
    Nmuabove100=np.zeros((len(Ei)-1,len(theta_earth),len(x)))
 
    for i in range(len(Ei)-1):
        for t in range(len(theta_earth)):
            j=0
            while (E[i][j]>1.2 and x[j]/2.65/100<14000):
                if (E[i][j]>1*10**2):
                    Nmuabove100[i][t][j]=Nmu[t][i]*2*np.pi
                j=j+1

    #Building the mass attenuation function for each rock material
    (funcSi, en_Si)=InterpMassAtt("Si")
    (funcAl, en_Al)=InterpMassAtt("Al")
    (funcO, en_O)=InterpMassAtt("O")
    (funcFe, en_Fe)=InterpMassAtt("Fe")
    (funcCa, en_Ca)=InterpMassAtt("Ca")
    (funcMg, en_Mg)=InterpMassAtt("Mg")
    (funcNi, en_Ni)=InterpMassAtt("Ni")
    (funcI, en_I)=InterpMassAtt("I")
    (funcNa, en_Na)=InterpMassAtt("Na")
    
    en1_mantle = GatherEnergies(en_Al, en_Si)
    en2_mantle = GatherEnergies(en1_mantle, en_O)  
    en3_mantle = GatherEnergies(en2_mantle, en_Fe)
    en4_mantle = GatherEnergies(en3_mantle, en_Ca)
    enfinal_mantle = GatherEnergies(en4_mantle, en_Mg)

    enfinal_core = GatherEnergies(en_Fe, en_Ni)
    enfinal_NaI = GatherEnergies(en_Na, en_I)

    att_mantle=np.zeros(len(enfinal_mantle))  
    for i in range(len(enfinal_mantle)):#doi:10.1016/j.epsl.2004.12.005
        att_mantle[i]= 0.2*funcSi(enfinal_mantle[i])+0.448*funcO(enfinal_mantle[i])+0.03*funcAl(enfinal_mantle[i])
        +0.08*funcFe(enfinal_mantle[i])+0.014*funcCa(enfinal_mantle[i])+0.228*funcMg(enfinal_mantle[i])
    funcmantle=inter.interp1d(enfinal_mantle, att_mantle*3.67*(1.973*10**(-14)), kind='linear', fill_value="extrapolate") #GeV, with energy given in MeV # Mass attenuatioin for the rock.

    att_core=np.zeros(len(enfinal_core))  
    for i in range(len(enfinal_core)):#IAG (Eder Molina:www.iag.usp.br/~eder/estrutura_e_evolucao.pdf)
        att_core[i]= 0.9*funcFe(enfinal_core[i])+0.1*funcNi(enfinal_core[i])
    funccore=inter.interp1d(enfinal_core, att_core*3.67*(1.973*10**(-14)), kind='linear', fill_value="extrapolate") #GeV, with energy given in MeV # Mass attenuatioin for the rock.

    att_NaI=np.zeros(len(enfinal_NaI))  
    for i in range(len(enfinal_NaI)):
        att_NaI[i]= 0.5*funcNa(enfinal_NaI[i])+0.5*funcI(enfinal_NaI[i])
    funcNaI=inter.interp1d(enfinal_NaI, att_NaI*3.67*(1.973*10**(-14)), kind='linear', fill_value="extrapolate") #GeV, with energy given in MeV # Mass attenuatioin for the rock.

    # Calculate how many DPs interact directly with the crystals.
    (NDPs_tot, Ndec_tot, Ndec_res, totalDPs, NOsc_tot) = CalcEvents(Ei, theta_earth, E, mDP, x, kappas, Nmuabove100, Na, epsilon, funcmantle, funccore, funcNaI, R_earth)

    N = len(NDPs_tot[0][0])
    # NCompsNaI = np.zeros(N)
    # NAbsNaI = np.zeros(N)
    # NOscNaI = np.eros(N)
    NdecNaI = np.zeros(N)
    NDPsNaI = np.zeros(N)
    Ndec = np.zeros(N)

    for i in range(0, N):
        # NCompsNaI[i] = round(np.sum(NComps_tot[:,:,i]),4)
        # NAbsNaI[i] = round(np.sum(NAbs_tot[:,:,i]),4)
        # NOscNaI[i] = round(np.sum(NOsc_tot[:,:,i]),4)
        NdecNaI[i] = round(np.sum(Ndec_tot[:,:,i]),4)
        NDPsNaI[i] = round(np.sum(NDPs_tot[:,:,i]),4)
        Ndec[i] = round(np.sum(Ndec_res[:,:,i]),4)

    if (mDP>1E-03):
        E_z = np.linspace(mDP, 3, N) #MeV
                    
    else:
        E_z = np.linspace(1E-03, 3, N) #MeV
      
    # with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/NDet_int_m=%.8f_k=%.11f.txt' %(mDP,kappas), 'w') as archive:
    #     for i in range(0, N):
    #         archive.write('{:.11f}  {:.16f}  {:.16f}  {:.16f}  {:.16f}\n'.format(E_z[i], NCompsNaI[i], NAbsNaI[i], NOscNaI[i], NdecNaI[i]))

    # with open('/home/dfreitas/DPnum_Comp/DirDetec/Spectra/NDir_int_m=%.8f_k=%.11f.txt' %(mDP,kappas), 'w') as archive:
    #     for i in range(0, N):
    #         archive.write('{:.11f}  {:.16f}  {:.16f}  {:.16f}  {:.16f}\n'.format(E_z[i], NCompsNaI[i], NAbsNaI[i], NOscNaI[i], NdecNaI[i]))

    print("total DPs = \n", np.sum(NDPs_tot))
    print("\nNumber of Bremsstrahlung producing DPs: \n", totalDPs)
    print("\nNumber of detections:\n", np.sum(Ndec_tot))
    print("\nNumber of events inside Earth:\n", np.sum(Ndec_res))
    print("\nNumber of oscillations inside crystals: \n", np.sum(NOsc_tot))
    #print("\nNumber of Comptons processes: \n", "{:e}".format(np.sum(NComps_res)))
    #print("\nPorcentagem de DP's que sofrem Compton: %d" %(np.sum(NComps_res)/np.sum(NDPs_tot)))

main()