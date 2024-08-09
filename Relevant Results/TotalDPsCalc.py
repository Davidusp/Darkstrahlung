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
import numba

epsilon=np.logspace(-10,-1,49)
mDP=np.array([0.000001,0.000002,0.000003,0.000004,0.000005,0.000007,0.000010,0.000013,0.000017,0.000022,0.000029,0.000039,
              0.000052,0.000069,0.000091,0.000121,0.000160,0.000212,0.000281,0.000373,0.000494,0.000655,0.000869,0.001151,
              0.001526,0.0011514,0.002024,0.002683,0.003556,0.004715,0.006251,0.008286,0.010985,0.014563,0.019307,0.025595,
              0.033932,0.044984,0.059636,0.079060,0.104811,0.138950,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.000000]) #MeV

totalDPs=np.zeros((len(mDP),len(epsilon)))

for z in range(len(mDP)):
    for k in range(len(epsilon)):
        
        arq1 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Absorption/Absorption Number/DPtot_Abs_m=%f00_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        texto1 = arq1.readlines() #salva cada linha do arquivo como uma string
        arq1.close()

        arq2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/DPs-Comp_PS/Comp Only/NComp_m=%f00_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        texto2 = arq2.readlines() #salva cada linha do arquivo como uma string
        arq2.close()

        arq3 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Nbremss/Results2/NBremss_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        texto3 = arq3.readlines() #salva cada linha do arquivo como uma string
        arq3.close()

        for i in range(len(texto1)): #separa as colunas e elimina a ultima coluna que contem '\n'
            texto1[i] = texto1[i].split("   ")
            texto1[i] = [float(x) for x in texto1[i]] #transforma os valores em numero
            texto2[i] = texto2[i].split("   ")
            texto2[i] = [float(x) for x in texto2[i]] #transforma os valores em numero
            texto3[i] = texto3[i].split("   ")
            texto3[i] = [float(x) for x in texto3[i]] #transforma os valores em numero 

        N=len(texto1)

        for i in range(N):
            totalDPs[z][k] = texto3[i][1] - texto1[i][1] - (13/64)*texto2[i][1]

ntotalDPstrans=np.transpose(totalDPs)

for i in range(len(mDP)):
    mDP[i]=1E+6*mDP[i]

# plt.pcolormesh(mDP,epsilon,ntotalDPstrans,cmap='viridis',norm=mc.LogNorm(), shading='nearest')
plt.pcolormesh(mDP,epsilon,ntotalDPstrans,cmap='viridis',norm=mc.LogNorm())
plt.colorbar()
plt.title('Total DPs -direct int ($year^{-1}$)',fontsize=20)
plt.xlabel('DP Mass (eV)',fontsize=14,fontweight='bold')
plt.ylabel('$\epsilon$',fontsize=18,fontweight='bold')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.xlim(1,1E+6)
plt.savefig('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/NtotDPs.png',dpi=100)