from cProfile import label
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import scipy.interpolate as inter


def main():

    epsilon=np.logspace(-10,-1,49)
    mDP=np.array([0.000001,0.000002,0.000003,0.000004,0.000005,0.000007,0.000010,0.000013,0.000017,0.000022,0.000029,0.000039,
              0.000052,0.000069,0.000091,0.000121,0.000160,0.000212,0.000281,0.000373,0.000494,0.000655,0.000869,0.001151,
              0.001526,0.0011514,0.002024,0.002683,0.003556,0.004715,0.006251,0.008286,0.010985,0.014563,0.019307,0.025595,
              0.033932,0.044984,0.059636,0.079060,0.104811,0.138950,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.000000]) #MeV 

    crystalMass = 106 #kg
    binsize = 2 #keV

    Ebins = np.linspace(0, 20000, 10001)
    energies = np.linspace(0,20000,60001) #keV

    for z in range(len(mDP)):
        for k in range(len(epsilon)):
            epsilon[k] = 0.01154781985
            mDP[z] = 0.75431200 
            arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/Final/Right Energy/NDir_int_m=%.6f00_k=%.11f.txt' %(mDP[z], epsilon[k]),'r')

            E_k = []
            Comp = []
            Abs = []
            Osc = []

            for linha in arq:
                line = linha.split("  ")
                E_k.append(float(line[0]))
                Comp.append(float(line[1]))
                Abs.append(float(line[2]))
                Osc.append(float(line[3]))
            arq.close()

            E_k = np.array(E_k)*1e03 #keV
            Comp = np.array(Comp)
            Abs = np.array(Abs)
            Osc = np.array(Osc)
    
            specComp = []
            specAbs = []
            specOsc = []

            Compfunc = inter.interp1d(E_k, Comp, kind = 'linear', fill_value="extrapolate")
            Absfunc = inter.interp1d(E_k, Abs, kind = 'linear', fill_value="extrapolate")
            Oscfunc = inter.interp1d(E_k, Osc, kind = 'linear', fill_value="extrapolate")

            sumComp = 0
            sumAbs = 0
            sumOsc = 0

            for i in range(1,len(energies)):
                if energies[i]%binsize != 0:
                    sumComp += Compfunc(energies[i])
                    sumAbs += Absfunc(energies[i])
                    sumOsc += Oscfunc(energies[i])

                else:
                    specComp.append(sumComp/5)
                    specAbs.append(sumAbs/5)
                    specOsc.append(sumOsc/5)
                    sumComp = 0
                    sumAbs = 0
                    sumOsc = 0

            specComp = np.array(specComp)/binsize/crystalMass/365.25
            specAbs = np.array(specAbs)/binsize/crystalMass/365.25
            specOsc = np.array(specOsc)/binsize/crystalMass/365.25

            #Plotting averages of each bin in the middle of the energy bin (1 keV by 1 keV)
            f1=plt.figure()
            plt.title("Dark photon spectrum - NaI(Tl). ($m_{Z'}=755 keV,\epsilon\sim 1.2\cdot 10^{-2}$)")
            # plt.plot(Ebins, specComp, c='blue', linewidth=2, label='Compton')
            # plt.plot(Ebins, specAbs, c='red',linestyle='--', linewidth=2, label='Absorption')
            plt.stairs(specComp, edges=Ebins, color='blue', label='Compton')
            plt.stairs(specAbs, edges=Ebins, color='red', label='Absorption')
            plt.stairs(specOsc, edges=Ebins, color='green', label='Oscillation in $\gamma$')
            plt.xlabel("$E_{Z'}[keV]$", fontsize=10)
            plt.ylabel("$Events~[counts/day/keV/kg]$", fontsize=12)
            plt.yscale('log')
            plt.xscale('linear')
            plt.xticks(fontsize=11)
            plt.yticks(fontsize=11)
            #plt.ylim(1E+2, 5e+9)
            #plt.xlim(1E-2,1)
            plt.legend(loc = 'lower center', fontsize=10)
            plt.grid(True)
            plt.savefig("/home/davidc/Documents/Master's Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/Final/Right Energy/Spec_int_m=%.8f_k=%.11f.png" %(mDP[z], epsilon[k]))
            #plt.show()
            plt.close()

            exit()

main()