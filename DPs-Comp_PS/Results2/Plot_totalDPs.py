import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from matplotlib.pyplot import figure
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.interpolate as inter
from scipy.integrate import quad
import matplotlib.colors as mc
import sys
#import numba
#from numba import jit

epsilon=np.logspace(-10,-1,49)
mDP=np.array([0.000001,0.000002,0.000003,0.000004,0.000005,0.000007,0.000010,0.000013,0.000017,0.000022,0.000029,0.000039,
              0.000052,0.000069,0.000091,0.000121,0.000160,0.000212,0.000281,0.000373,0.000494,0.000655,0.000869,0.001151,
              0.001526,0.0011514,0.002024,0.002683,0.003556,0.004715,0.006251,0.008286,0.010985,0.014563,0.019307,0.025595,
              0.033932,0.044984,0.059636,0.079060,0.104811,0.138950,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.000000]) #MeV

totalDPs=np.zeros((len(mDP),len(epsilon)))

for z in range(len(mDP)):
    for k in range(len(epsilon)):

        arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/DPs-Comp_PS/Results2/DPtot_Comp_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        texto = arq.readlines() #salva cada linha do arquivo como uma string
        arq.close()

        for i in range(len(texto)): #separa as colunas e elimina a ultima coluna que contem '\n'
            texto[i] = texto[i].split("   ")
            texto[i] = [float(x) for x in texto[i]] #transforma os valores em numero 

        N=len(texto)

        for i in range(N):
            totalDPs[z][k]=texto[i][1]

ntotalDPstrans=np.transpose(totalDPs)

for i in range(len(mDP)):
    mDP[i]=1E+6*mDP[i]

# plt.pcolormesh(mDP,epsilon,ntotalDPstrans,cmap='viridis',norm=mc.LogNorm(), shading='nearest')
plt.pcolormesh(mDP,epsilon,ntotalDPstrans,cmap='viridis',norm=mc.LogNorm())
plt.colorbar()
plt.title(' Total DPs (-Comp) ($year^{-1}$)',fontsize=20)
plt.xlabel('DP Mass (eV)',fontsize=14,fontweight='bold')
plt.ylabel('$\epsilon$',fontsize=18,fontweight='bold')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.xlim(1,1E+6)
plt.savefig('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/DPs-Comp_PS/Results2/NDPstot.png',dpi=100)