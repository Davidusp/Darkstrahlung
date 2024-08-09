import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


def main():

    arq = open('PeXS.txt','r')
    E_k = []
    sig = []

    for linha in arq:
        line = linha.split()
        E_k.append(float(line[0]))
        sig.append(float(line[1]))
    arq.close()

    print(E_k)
    
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
    #plt.savefig("/home/davidc/Documents/Master's Analisys/Results Compton/CS_%dE_%dcos_k=%f_mphi=%f.png" %(len(E_k), 15000,1e-4, i))
    plt.show()
    #plt.close()

main()