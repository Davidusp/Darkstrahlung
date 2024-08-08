from doctest import REPORT_CDIFF
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
#from numba import jit

# epsilon=np.array([0.0001,0.1,0.00015399265,0.00023713737,0.00036517413,0.00056234133,0.00086596432,0.00133352143,0.0177827941,
#                   0.00205352503,0.00316227766,0.00486967525,0.00749894209,0.01154781985,0.02738419634,0.06493816316,1.5e-10,
#                   1.33e-09,1.155e-08,1.778e-08,1.5399e-07,1.33352e-06,1.154782e-05,1.778279e-05,1e-07,1e-10,2.4e-10,2.05e-09,
#                   2.738e-08,2.3714e-07,2.05353e-06,2.73842e-05,3.7e-10,3.16e-09,3.6517e-07,3.16228e-06,4.87e-09,4.217e-08,
#                   4.86968e-06,4.216965e-05,5.6e-10,5.6234e-07,6.494e-08,6.493816e-05,7.5e-09,7.49894e-06,8.7e-10,8.6596e-07])

# mDP=np.array([0.00016,0.000121,0.000212,0.000281,0.000373,0.000494,0.000655,0.000869,0.001151,0.001526,0.002024,0.002683,0.003556,
#               0.004715,0.006251,0.07906,0.008286,0.010985,0.0011514,0.13895,0.014563,0.019307,0.025595,0.033932,0.044984,0.059636,
#               0.104811,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.0,1.3e-05,1.7e-05,1e-05,1e-06,2.2e-05,2.9e-05,2e-06,
#               3.9e-05,3e-06,4e-06,5.2e-05,5e-06,6.9e-05,7e-06,9.1e-05]) #MeV

#epsilon=np.logspace(-10,-1,49)
mDP=np.array([0.000001,0.000002,0.000003,0.000004,0.000005,0.000007,0.000010,0.000013,0.000017,0.000022,0.000029,0.000039,
              0.000052,0.000069,0.000091,0.000121,0.000160,0.000212,0.000281,0.000373,0.000494,0.000655,0.000869,0.001151,
              0.001526,0.0011514,0.002024,0.002683,0.003556,0.004715,0.006251,0.008286,0.010985,0.014563,0.019307,0.025595,
              0.033932,0.044984,0.059636,0.079060,0.104811,0.138950,0.184207,0.244205,0.323746,0.429193,0.568987,0.754312,1.000000]) #MeV

epsilon = np.array([0.00000000010,0.00000000015,0.00000000024,0.00000000037,0.00000000056,0.00000000087,0.00000000133,0.00000000205,
           0.00000000316,0.00000000487,0.00000000750,0.00000001155,0.00000001778,0.00000002738,0.00000004217,0.00000006494,
           0.00000010000,0.00000015399,0.00000023714,0.00000036517,0.00000056234,0.00000086596,0.00000133352,0.00000205353,
           0.00000316228,0.00000486968,0.00000749894,0.00001154782,0.00001778279,0.00002738420,0.00004216965,0.00006493816,
           0.00010000000,0.00015399265,0.00023713737,0.00036517413,0.00056234133,0.00086596432,0.00133352143,0.00205352503,
           0.00316227766,0.00486967525,0.00749894209,0.01154781985,0.01778279410,0.02738419634,0.04216965034,0.06493816316,
           0.10000000000])

# mDP = np.sort(mDP)
# epsilon = np.sort(epsilon)

totalOsc=np.zeros((len(mDP),len(epsilon)))
totalBremss = np.zeros((len(mDP),len(epsilon)))

for z in range(len(mDP)):
    for k in range(len(epsilon)):
        
        #print(epsilon[k])
        arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Oscillation/Results2/Final/noscRock_m=%f00_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        arq2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Nbremss/Results2/NBremss_m=%f00_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
        texto = arq.readlines() #salva cada linha do arquivo como uma string
        texto2 = arq2.readlines()
        arq.close()
        arq2.close()

        for i in range(len(texto)): #separa as colunas e elimina a ultima coluna que contem '\n'
            texto[i] = texto[i].split("   ")
            texto2[i] = texto2[i].split("   ")
            texto[i] = [float(x) for x in texto[i]] #transforma os valores em numero 
            texto2[i] = [float(x) for x in texto2[i]]

        N=len(texto)

        for i in range(N):
            totalOsc[z][k]=texto[i][1]
            totalBremss[z][k]=texto2[i][1]
            # if totalOsc[z][k] != 0:
            #     print(totalOsc[z][k])

ntotalDPstrans=np.transpose(totalOsc)
nbremsstrans=np.transpose(totalBremss)

for i in range(len(mDP)):
    mDP[i]=1E+6*mDP[i]

#plt.pcolormesh(mDP,epsilon,ntotalDPstrans,cmap='viridis',norm=mc.LogNorm(), shading='nearest')
plt.pcolormesh(mDP,epsilon,ntotalDPstrans/nbremsstrans*100,cmap='viridis',norm=mc.LogNorm())
plt.colorbar()
plt.title('Oscillations/total DPs ($year^{-1}$)',fontsize=20)
plt.xlabel('DP Mass (eV)',fontsize=14,fontweight='bold')
plt.ylabel('$\epsilon$',fontsize=18,fontweight='bold')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
# plt.xlim(1,1E+6)
# plt.ylim(1.7*1e-7,1e-1)
plt.savefig('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Oscillation/OscNumber.png',dpi=100)