from cProfile import label
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import scipy.interpolate as inter
import pandas as pd


def main():

    epsilon=np.logspace(-10,-1,49)
    mDP=np.array([0.104811,0.138950,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.000000,1.2,
                  1.57079914,2.05617494,2.69153152,3.52321282,4.61188305,6.03695159,7.90236529,10.34419048,
                  13.54053789,17.72455436,23.20142891,30.37065375,39.75516391,52.03948096,68.11964314,
                  89.16856387,116.72158596,152.78847205,200.]) #MeV
    
    # mDP=np.array([0.104811,0.138950,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.000000,1.2,1.57079914,
    #               2.05617494,2.69153152,3.52321282,4.61188305,6.03695159,7.90236529,10.34419048,13.54053789,
    #               17.72455436]) #MeV
    
    crystalMass = 106 #kg
    binsize = 2 #keV
    #PLOTAR ATÉ 5*4 MeV -> máxima resolução dos cristais p/ 5 cristais (SET3)

    Ebins = np.linspace(0, 100, 51)
    energies = np.linspace(0,100,301) #keV
    N=50

    for z in range(len(mDP)):
        for k in range(len(epsilon)):
            mDP = 0.00828600 
            epsilon = 0.00004216965
            arq2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/SpecAllEarth/100keV/NDir100keV_int_m=%.8f_k=%.11f.txt' %(mDP, epsilon),'r')

            E_k = []
            Ndet = []
            Ndec = []

            for linha in arq2:
                line = linha.split("  ")
                E_k.append(float(line[0]))
                
                if mDP <=1:
                    Ndet.append(float(line[1])+float(line[2]))
                    Ndec.append(float(line[3]))
                else:
                    Ndet.append(float(line[1]))
                    Ndec.append(float(line[2]))
                
            arq2.close()

            E_k = np.array(E_k)*1e03 #keV
            Ndet = np.array(Ndet)
            Ndec = np.array(Ndec)

            if mDP <=1:
                arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/Final/100keV/NDir100keV_int_m=%.8f_k=%.11f.txt' %(mDP, epsilon),'r')
                texto = arq.readlines()
                arq.close()

                for i in range(len(texto)):
                    texto[i] = texto[i].split("   ")
                    texto[i] = texto[i][0].split("  ")
                    for j in range(len(texto[i])):
                        texto[i][j] = float(texto[i][j])
                texto = np.array(texto)

                for l in range(len(Ndet)):
                    Ndet[l] += texto[l][1]+texto[l][2]+texto[l][3]

                NDPs2 = np.loadtxt('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/SpecAllEarth/DPsOnly/NDPstot_m=%.8f_k=%.11f.txt' %(mDP, epsilon))
                NDPs1 = np.loadtxt('/home/davidc/Documents/Master\'s Analisys/Parameter space/Nbremss/NBremss+energy/NBremss_m=%.13f_k=%.11f.txt' %(mDP, epsilon))
                NDPs1 = NDPs1[:,1]
                #NDPs2 = NDPs2[:,1]

                NDPs = NDPs1+NDPs2

                #print(E_k, Ndet)
            
            else:
                NDPs = np.loadtxt('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/SpecAllEarth/DPsOnly/NDPstot_m=%.8f_k=%.11f.txt' %(mDP, epsilon))

            specdet = []
            specdec = []
            specDPs = []

            #print(E_k, Ndet)

            Detfunc = inter.interp1d(E_k, Ndet, kind = 'linear', fill_value="extrapolate")
            Decfunc = inter.interp1d(E_k, Ndec, kind = 'linear', fill_value="extrapolate")
            #DPsfunc = inter.interp1d(E_k, NDPs, kind = 'linear', fill_value="extrapolate")

            sumDet = 0
            sumDec = 0
            #sumOsc = 0
            sumDPs = 0

            for i in range(1,len(energies)):
                if energies[i]%binsize != 0:
                    sumDet += Detfunc(energies[i])
                    sumDec += Decfunc(energies[i])
                    #sumDPs += DPsfunc(energies[i])

                else:
                    specdet.append(sumDet/5) #Computing the mean values in each energy bin
                    specdec.append(sumDec/5)
                    #specOsc.append(sumOsc/5)
                    specDPs.append(sumDPs/5)
                    sumDet = 0
                    sumDec = 0
                    #sumOsc = 0
                    sumDPs = 0

            specdet = np.array(specdet)/binsize/crystalMass/365.25
            specdec = np.array(specdec)/binsize/crystalMass/365.25
            #specOsc = np.array(specOsc)/binsize/crystalMass/365.25
            #specDPs = np.array(specDPs)/binsize/crystalMass/365.25

            x_start = 0
            x_end = mDP*1000

            #Plotting averages of each bin in the middle of the energy bin (1 keV by 1 keV)
            f1=plt.figure()
            plt.title("Dark photon spectrum ($m_{Z'}=8.3 keV,\epsilon\sim 4.2\cdot 10^{-5}$)")
            # plt.plot(Ebins, specComp, c='blue', linewidth=2, label='Compton')
            # plt.plot(Ebins, specAbs, c='red',linestyle='--', linewidth=2, label='Absorption')
            plt.stairs(specdet, edges=Ebins, color='blue', label='Detections (NaI)')
            plt.stairs(specdec, edges=Ebins, color='red', label='Interactions with Earth')
            #plt.stairs(specOsc, edges=Ebins, color='green', label='Oscillation in $\gamma$')
            #plt.stairs(specDPs, edges=Ebins, color='yellow', label='total DPs')

            # Plot the shaded area
            plt.axvspan(x_start, x_end, color='gray', alpha=0.3)

            plt.xlabel("$E_{Z'}[keV]$", fontsize=11)
            plt.ylabel("$Events~[counts/day/keV/kg]$", fontsize=12)
            plt.yscale('log')
            plt.xscale('linear')
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=11)
            plt.legend(loc = 'lower right', fontsize="9")
            plt.grid(True)
            plt.savefig("/home/davidc/Documents/Master's Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/SpecAllEarth/DetSpectra100keV/NDet100keV_m=%.8f_k=%.11f.png" %(mDP, epsilon))
            #plt.show()
            plt.close()

            exit()

main()