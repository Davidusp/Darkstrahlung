import numpy as np # type: ignore
import matplotlib # type: ignore
import matplotlib.pyplot as plt # type: ignore
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
mDP=np.array([0.000001,0.000002,0.000003,0.000004,0.000005,0.000007,0.000010,0.000013,0.000017,0.000022,0.000029,
              0.000039,0.000052,0.000069,0.000091,0.000121,0.000160,0.000212,0.000281,0.000373,0.000494,0.000655,
              0.000869,0.001151,0.001526,0.0011514,0.002024,0.002683,0.003556,0.004715,0.006251,0.008286,0.010985,
              0.014563,0.019307,0.025595,0.033932,0.044984,0.059636,0.079060,0.104811,0.138950,0.184207,0.244205,
              0.323746,0.429193,0.568987,0.754312,1.000000,1.2,1.57079914,2.05617494,2.69153152,3.52321282,4.61188305,
              6.03695159,7.90236529,10.34419048,13.54053789,17.72455436,23.20142891,30.37065375,39.75516391,52.03948096,
              68.11964314,89.16856387,116.72158596,152.78847205,200.]) #MeV

totalDet = np.zeros((len(mDP),len(epsilon)))
totalBremss = np.zeros((len(mDP),len(epsilon)))

for z in range(len(mDP)):
    for k in range(len(epsilon)):
        
        if mDP[z] <= 1:
            arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/Final/Right Energy/NDir_int_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
            arqbremss1 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Nbremss/Results2/NBremss_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
            arq2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/SpecAllEarth/NDir_int_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r') #abre o arquivo
            arqbremss2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/SpecAllEarth/DPsOnly/NDPstot_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r')

            texto = arq.readlines() #salva cada linha do arquivo como uma string
            textobremss1 = arqbremss1.readlines()
            texto2 = arq2.readlines()
            textobremss2 = arqbremss2.readlines()

            arq.close()
            arqbremss1.close()
            arq2.close()
            arqbremss2.close()

            textobremss1[0] = textobremss1[0].split("   ")
            textobremss1[0] = [float(x) for x in textobremss1[0]]
            textobremss1 = textobremss1[0]

            for i in range(len(texto)): #separa as colunas e elimina a ultima coluna que contem '\n'
                texto[i] = texto[i].split("   ")
                texto[i] = texto[i][0].split("  ")
                texto2[i] = texto2[i].split("   ")
                texto2[i] = texto2[i][0].split("  ")
                textobremss2[i] = float(textobremss2[i])
                for j in range(len(texto[i])):
                    texto[i][j] = float(texto[i][j])
                    texto2[i][j] = float(texto2[i][j])

            texto = np.array(texto)
            texto2 = np.array(texto2)
            textobremss1 = np.array(textobremss1)

            totalDet[z][k]=np.sum(texto[:,1])+np.sum(texto[:,2])+np.sum(texto[:,3])+np.sum(texto2[:,1])+np.sum(texto2[:,2])
            totalBremss[z][k]=textobremss1[1]+np.sum(textobremss2)
            
        else:
            arq2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/SpecAllEarth/NDir_int_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r')
            arqbremss2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/SpecAllEarth/DPsOnly/NDPstot_m=%.8f_k=%.11f.txt' %(mDP[z], epsilon[k]), 'r')
        
            texto2 = arq2.readlines()
            textobremss2 = arqbremss2.readlines()
        
            arq2.close()
            arqbremss2.close()

            for i in range(len(texto2)): #separa as colunas e elimina a ultima coluna que contem '\n'
                texto2[i] = texto2[i].split("   ")
                texto2[i] = texto2[i][0].split("  ")
                textobremss2[i] = float(textobremss2[i])
                for j in range(len(texto2[i])):
                    texto2[i][j] = float(texto2[i][j])
            
            texto2 = np.array(texto2)
            textobremss2 = np.array(textobremss2)

            totalDet[z][k]=np.sum(texto2[:,1])
            #if totalDet[z][k] > 0:
                #print(totalDet[z][k])
            totalBremss[z][k]=np.sum(textobremss2)

fractionNum = np.zeros((len(mDP),len(epsilon)))

for z in range(len(totalDet)):
    for k in range(len(totalDet[z])):
        if totalDet[z][k] < 0:
            totalDet[z][k] = 0
        if totalBremss[z][k] != 0:
            fractionNum[z][k] = totalDet[z][k]/totalBremss[z][k]

ntransDetect = np.transpose(fractionNum)

limitmDP = []
limitepsilon = []

for z in range(len(totalDet)):
    for k in range(len(totalDet[z])):
        if int(totalDet[z][k]) == 1:
            limitmDP.append(mDP[z]*1e06)
            limitepsilon.append(epsilon[k])

#funclimit = inter.interp1d(limitmDP, limitepsilon, kind='linear',fill_value="extrapolate")

for i in range(len(mDP)):
    mDP[i]=1E+6*mDP[i]

# Create a figure and axes
fig, ax = plt.subplots()

# Plot 2D plot with pcolormesh
pcm = ax.pcolormesh(mDP, epsilon, ntransDetect, cmap='viridis', norm=mc.LogNorm())

# Plot scatter plot
ax.scatter(limitmDP, limitepsilon, c='red', s=4)

# Add colorbar
plt.colorbar(pcm, ax=ax)

# Set labels and title
ax.set_title('Detections/Bremsstrahlung ($year^{-1}$)', fontsize=14)
ax.set_xlabel('DP Mass (eV)', fontsize=14, fontweight='bold')
ax.set_ylabel('$\epsilon$', fontsize=18, fontweight='bold')
ax.tick_params(axis='both', which='minor', labelsize=16)
# ax.xaxis.set_tick_params(labelsize=12.5)
# ax.yaxis.set_tick_params(labelsize=12.5)
ax.set_xscale('log')
ax.set_yscale('log')

# Save the figure
plt.savefig('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/NDetection.png', dpi=100)
#plt.show()


# # plt.pcolormesh(mDP,epsilon,ntotalDPstrans,cmap='viridis',norm=mc.LogNorm(), shading='nearest')
# plt.pcolormesh(mDP,epsilon,ntransDetect,cmap='viridis',norm=mc.LogNorm())
# plt.scatter(limitmDP,limitepsilon,c='red',s=4)
# plt.colorbar()
# plt.title('Detections ($year^{-1}$)',fontsize=14)
# plt.xlabel('DP Mass (eV)',fontsize=14,fontweight='bold')
# plt.ylabel('$\epsilon$',fontsize=18,fontweight='bold')
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.xscale('log')
# plt.yscale('log')
# #plt.xlim(1,1E+6)
# #plt.show()
# plt.savefig('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/NDetectionTot.png',dpi=100)