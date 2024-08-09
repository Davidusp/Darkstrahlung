from cProfile import label
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import scipy.interpolate as inter


def main():

    arq = open('NDir_int_m=0.00012100_k=0.00133352143.txt','r')
    E_k = []
    Comp = []
    Abs = []

    for linha in arq:
        line = linha.split("  ")
        E_k.append(float(line[0]))
        Comp.append(float(line[1]))
        Abs.append(float(line[2]))
    arq.close()

    E_k = np.array(E_k)*1e03 #keV
    Comp = np.array(Comp)
    Abs = np.array(Abs)

    crystalMass = 106 #kg
    binsize = 2 #keV

    #Ebins = np.linspace(1,2999,1500)
    Ebins = np.linspace(0, 3000, 1501)
    energies = np.linspace(0,3000,9001)
    specComp = []
    specAbs = []

    Compfunc = inter.interp1d(E_k, Comp, kind = 'linear', fill_value="extrapolate")
    Absfunc = inter.interp1d(E_k, Abs, kind = 'linear', fill_value="extrapolate")

    sumComp = 0
    sumAbs = 0

    for i in range(1,len(energies)):
        if energies[i]%binsize != 0:
            sumComp += Compfunc(energies[i])
            sumAbs += Absfunc(energies[i])

        else:
            specComp.append(sumComp/4)
            specAbs.append(sumAbs/4)
            sumComp = 0
            sumAbs = 0

    specComp = np.array(specComp)/binsize/crystalMass/365.25
    specAbs = np.array(specAbs)/binsize/crystalMass/365.25

    #Plotting averages of each bin in the middle of the energy bin (1 keV by 1 keV)
    f1=plt.figure()
    plt.title('Direct detection in NaI(Tl) crystals')
    # plt.plot(Ebins, specComp, c='blue', linewidth=2, label='Compton')
    # plt.plot(Ebins, specAbs, c='red',linestyle='--', linewidth=2, label='Absorption')
    plt.stairs(specComp, edges=Ebins, color='blue', label='Compton')
    plt.stairs(specAbs, edges=Ebins, color='red', label='Absorption')
    plt.xlabel("$E_{k}[keV]$")
    plt.ylabel("$counts~[counts/day/keV/kg]$")
    plt.yscale('log')
    plt.xscale('linear')
    #plt.ylim(1E+2, 5e+9)
    #plt.xlim(1E-2,1)
    plt.legend(loc = 'upper right')
    plt.grid(True)
    #plt.savefig("/home/davidc/Documents/Master's Analisys/Results Compton/CS_%dE_%dcos_k=%f_mphi=%f.png" %(len(E_k), 15000,1e-4, i))
    plt.show()
    #plt.close()

main()