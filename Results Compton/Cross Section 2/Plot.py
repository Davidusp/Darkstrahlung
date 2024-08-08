import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


def main():

    m = [0.0000010,0.0000013,0.0000018,0.0000023,0.000003,0.000004,0.000005,0.000007,0.000010,0.000013,0.000017,0.000022,
         0.000029,0.000039,0.000052,0.000069,0.000091,0.000121,0.000160,0.000212,0.000281,0.000373,0.000494,0.000655,0.000869,
         0.001151,0.001526,0.002024,0.002683,0.003556,0.004715,0.006251,0.008286,0.010985,0.014563,0.019307,0.025595,0.033932,
         0.044984,0.059636,0.079060,0.104811,0.138950,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.000000]

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
        plt.title('Compton-like Cross-section')
        plt.plot(E_k, sig, c='blue', linewidth=2)
        plt.xlabel("$E_{k}[MeV]$", fontsize=13)
        plt.ylabel("$\sigma~[barn/atom]$",fontsize=13)
        plt.yscale('log')
        plt.xscale('linear')
        plt.xticks(fontsize=13)
        plt.yticks(fontsize=13)
        #plt.ylim(1E+2, 5e+9)
        #plt.xlim(1E-2,1)
        #plt.legend(loc = 'upper right')
        plt.grid(True)
        plt.savefig("/home/davidc/Documents/Master's Analisys/Results Compton/Cross Section 2/CS_%dE_%dcos_k=%f_mphi=%f.png" %(len(E_k), 15000,1e-4, i))
        #plt.show()
        plt.close()

main()