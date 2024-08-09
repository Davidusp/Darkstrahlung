import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('Agg')
from matplotlib.pyplot import figure
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.interpolate as inter
from scipy.integrate import quad
import matplotlib.colors as mc
import sys
import numba
#from numba import jit

epsilon=np.logspace(-10,-1,49)
mDP=np.array([0.000001,0.000002,0.000003,0.000004,0.000005,0.000007,0.000010,0.000013,0.000017,0.000022,0.000029,0.000039,
              0.000052,0.000069,0.000091,0.000121,0.000160,0.000212,0.000281,0.000373,0.000494,0.000655,0.000869,0.001151,
              0.001526,0.0011514,0.002024,0.002683,0.003556,0.004715,0.006251,0.008286,0.010985,0.014563,0.019307,0.025595,
              0.033932,0.044984,0.059636,0.079060,0.104811,0.138950,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.000000]) #MeV

# mDP=np.array([0.0000000001000,0.000000001668,0.000000002783,0.000000004642,0.000000007743,0.000000012915,0.000000021544,
#               0.000000035938,0.000000059948,0.000000100000]) #MeV

totalCmp = np.zeros((len(mDP),len(epsilon)))
totalAbs = np.zeros((len(mDP),len(epsilon)))
totalBremss = np.zeros((len(mDP),len(epsilon)))
totalOsc = np.zeros((len(mDP),len(epsilon)))

for z in range(len(mDP)):
    for k in range(len(epsilon)):

        # arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/DPs-Comp_PS/Comp Only/Results+/NComp_m=%f00_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        # arq2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Nbremss/Results2/NBremss_m=%f00_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        # arq3 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/Absorption Number/Results2/DPtot_Abs_m=%f00_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r')
        # arq4 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Oscillation/Results2/noscRock_m=%f00_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r')
        arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/DPs-Comp_PS/Comp Only/Results+/Final/NComp_m=%.13f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        arq2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Nbremss/Results2/NBremss_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        arq3 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/Absorption Number/Final/DPtot_Abs_m=%.13f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r')
        arq4 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Oscillation/Results2/Final/noscRock_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r')

        texto = arq.readlines() #salva cada linha do arquivo como uma string
        texto2 = arq2.readlines()
        texto3 = arq3.readlines()
        texto4 = arq4.readlines()
        arq.close()
        arq2.close()
        arq3.close()
        arq4.close()

        for i in range(len(texto)): #separa as colunas e elimina a ultima coluna que contem '\n'
            texto[i] = texto[i].split("   ")
            texto2[i] = texto2[i].split("   ")
            texto3[i] = texto3[i].split("   ")
            texto4[i] = texto4[i].split("   ")
            texto[i] = [float(x) for x in texto[i]] #transforma os valores em numero
            texto2[i] = [float(x) for x in texto2[i]]
            texto3[i] = [float(x) for x in texto3[i]]
            texto4[i] = [float(x) for x in texto4[i]]

        N=len(texto)

        for i in range(N):
            totalCmp[z][k]=texto[i][1]
            totalBremss[z][k]=texto2[i][1]
            totalAbs[z][k]=texto3[i][1]
            totalOsc[z][k]=texto4[i][1]
            # print("Nbremss = %s    NAbs = %s    NComp = %s    NOsc = %s\n" %("{:.20e}".format(totalBremss[z][k]),
            #                                                                  "{:.20e}".format(totalAbs[z][k]),
            #                                                                  "{:.20e}".format(totalCmp[z][k]),
            #                                                                  "{:.20e}".format(totalOsc[z][k])))
            
nDetect = totalBremss-totalAbs-totalCmp-totalOsc

for z in range(len(nDetect)):
    for k in range(len(nDetect[z])):
        if nDetect[z][k] < 0:
            nDetect[z][k] = 0

ntransDetect = np.transpose(nDetect)

limitmDP = []
limitepsilon = []

for z in range(len(totalBremss)):
    for k in range(len(totalBremss[z])):
        if int(totalBremss[z][k]) == 1:
            limitmDP.append(mDP[z]*1e06)
            limitepsilon.append(epsilon[k])

funclimit = inter.interp1d(limitmDP, limitepsilon, kind='linear',fill_value="extrapolate")

for i in range(len(mDP)):
    mDP[i]=1E+6*mDP[i]

#print(len(limitmDP), '\n', len(limitepsilon))

# plt.pcolormesh(mDP,epsilon,ntotalDPstrans,cmap='viridis',norm=mc.LogNorm(), shading='nearest')
plt.pcolormesh(mDP,epsilon,ntransDetect,cmap='viridis',norm=mc.LogNorm())
plt.plot(mDP,funclimit(mDP),c='red',linestyle='-')
plt.colorbar()
plt.title('Total DPs - extinction ($year^{-1}$)',fontsize=20)
plt.xlabel('DP Mass (eV)',fontsize=14,fontweight='bold')
plt.ylabel('$\epsilon$',fontsize=18,fontweight='bold')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
#plt.xlim(1,1E+6)
#plt.show()
plt.savefig('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Ntot-int.png',dpi=100)