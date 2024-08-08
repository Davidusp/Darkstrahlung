import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


def main():

    arq = open('PeXS_E vs Sig.txt','r') #Energy (MeV) x CS (cm²/g)
    E_k = []
    sig = []

    mDP = 0.000001 #MeV
    kappas = 0.06493816316
    q = 0.30282212088 #Elementary charge in natural units
    alpha = q**2/(4*np.pi) #Fine structure constant
    alphaline = (q*kappas)**2/(4*np.pi)
    for linha in arq:
        line = linha.split()
        E_k.append(float(line[0]))
        sig.append(float(line[1]))
    arq.close()

    for l in range(0, len(E_k)):
        sig[l] = (sig[l]/np.sqrt(1-mDP**2/E_k[l]**2))*(alphaline/alpha)*(1e24/2.74e22)

    sig = np.array(sig)

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
    plt.close()

main()