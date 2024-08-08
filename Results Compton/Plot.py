import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


def main():

    m = [0.000001,0.00000133,0.00000176,0.00000233,0.00000309,0.00000409,0.00000543,0.0000072,0.00000954,0.00001265,0.00001677,0.00002223,0.00002947,0.00003907,0.00005179,0.00006866,0.00009103,0.00012068,0.00015999,0.0002121,0.00028118,0.00037276,0.00049417,0.00065513,0.00086851,0.0011514,0.00152642,0.00202359,0.0026827,0.00355648,0.00471487,0.00625055,0.00828643,0.01098541,0.01456348,0.01930698,0.02559548,0.03393222,0.04498433,0.05963623,0.07906043,0.10481131,0.13894955,0.184207,0.24420531,0.32374575,0.42919343,0.5689866,0.75431201,1.]

    for i in m:
        arq = open('CSxE_m=%f_k=0.000100.txt' %(i),'r')
        E_k = []
        sig = []

        for linha in arq:
            line = linha.split()
            E_k.append(float(line[0]))
            sig.append(float(line[1]))
        arq.close()

        f1=plt.figure()
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
        plt.savefig("/home/davidc/Documents/Master's Analisys/Results Compton/CS_%dE_%dcos_k=%f_mphi=%f.png" %(len(E_k), 15000,1e-4, i))
        #plt.show()
        plt.close()

main()