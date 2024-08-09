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
import numba
#from numba import jit

epsilon=np.logspace(-4,-2,20)
mDP=np.logspace(-10,-3,30) #GeV

nbremss=np.zeros((len(mDP),len(epsilon)))

for z in range(len(mDP)):
    for k in range(len(epsilon)):
        arq = open('/home/LuisFranca/forDavid/times/nbremss/nbremss_%s_%s.txt' % ("{:.3e}".format(mDP[z]),"{:.3e}".format(epsilon[k])), 'r') #abre o arquivo
        texto = arq.readlines() #salva cada linha do arquivo como uma string
        arq.close()

        for i in range(len(texto)): #separa as colunas e elimina a ultima coluna que contem '\n'
           texto[i] = texto[i].split(" ")
           texto[i] = [float(x) for x in texto[i]] #transforma os valores em numero 

        N=len(texto)

        for i in range(N):
            nbremss[z][k]=texto[i][0]
                
nbremsstrans=np.transpose(nbremss)
for i in range(len(mDP)):
    mDP[i]=1E+9*mDP[i]
    
plt.pcolormesh(mDP,epsilon,nbremsstrans,cmap='viridis',norm=mc.LogNorm())
plt.colorbar()
plt.title(' ($year^{-1}$)',fontsize=20)
plt.xlabel('DP Mass (eV)',fontsize=14,fontweight='bold')
plt.ylabel('$\epsilon$',fontsize=18,fontweight='bold')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.1,1E+6)
plt.savefig('nbremss.png',dpi=100)