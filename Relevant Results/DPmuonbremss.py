import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.interpolate as inter
import sys

#mDP=float(sys.argv[1])
#print(mDP)
mDP = 0.000323746 #GeV

arq = open('/home/LuisFranca/forDavid/muonflux.tsv', 'r')
texto1 = arq.readlines()
arq.close()

for i in range(len(texto1)):
    texto1[i] = texto1[i].split("\t")
    texto1[i] = [float(x) for x in texto1[i]]

N1=len(texto1)
p=np.zeros(N1)
y=np.zeros(N1)
E=np.zeros(N1)
mu=0.105658 #GeV
Na=6.02*10**23

for i in range(N1):
    p[i]=texto1[i][0]
    y[i]=texto1[i][1]
    E[i]=(p[i]**2+mu**2)**0.5
    
func=inter.interp1d(p,y/p**3)

arq = open('/home/LuisFranca/forDavid/muon_stoppower.tsv', 'r')
texto1 = arq.readlines()
arq.close()

for i in range(len(texto1)):
    texto1[i] = texto1[i].split("\t")
    texto1[i] = [float(x) for x in texto1[i]]

N1=len(texto1)
dEdx=np.zeros(N1)
Em=np.zeros(N1)

for i in range(N1):
    Em[i]=texto1[i][0]
    dEdx[i]=texto1[i][1]

funcdEdx=inter.interp1d(Em,dEdx)

arq = open('/home/LuisFranca/forDavid/muonangle_1GeV.tsv', 'r')
texto1 = arq.readlines()
arq.close()

for i in range(len(texto1)):
    texto1[i] = texto1[i].split("\t")
    texto1[i] = [float(x) for x in texto1[i]]

N1=len(texto1)
zenite_1GeV=np.zeros(N1)
ratio_1GeV=np.zeros(N1)

for i in range(N1):
    zenite_1GeV[i]=np.arccos(texto1[i][0])
    ratio_1GeV[i]=texto1[i][1]
func1GeV=inter.interp1d(zenite_1GeV,ratio_1GeV)

arq = open('/home/LuisFranca/forDavid/muonangle_10GeV.tsv', 'r')
texto1 = arq.readlines()
arq.close()

for i in range(len(texto1)):
    texto1[i] = texto1[i].split("\t")
    texto1[i] = [float(x) for x in texto1[i]]

N1=len(texto1)
zenite_10GeV=np.zeros(N1)
ratio_10GeV=np.zeros(N1)

for i in range(N1):
    zenite_10GeV[i]=np.arccos(texto1[i][0])
    ratio_10GeV[i]=texto1[i][1]
func10GeV=inter.interp1d(zenite_10GeV,ratio_10GeV)

arq = open('/home/LuisFranca/forDavid/muonangle_100GeV.tsv', 'r')
texto1 = arq.readlines()
arq.close()

for i in range(len(texto1)):
    texto1[i] = texto1[i].split("\t")
    texto1[i] = [float(x) for x in texto1[i]]

N1=len(texto1)
zenite_100GeV=np.zeros(N1)
ratio_100GeV=np.zeros(N1)

for i in range(N1):
    zenite_100GeV[i]=np.arccos(texto1[i][0])
    ratio_100GeV[i]=texto1[i][1]
func100GeV=inter.interp1d(zenite_100GeV,ratio_100GeV)

arq = open('/home/LuisFranca/forDavid/muonangle_1000GeV.tsv', 'r')
texto1 = arq.readlines()
arq.close()

for i in range(len(texto1)):
    texto1[i] = texto1[i].split("\t")
    texto1[i] = [float(x) for x in texto1[i]]

N1=len(texto1)
zenite_1000GeV=np.zeros(N1)
ratio_1000GeV=np.zeros(N1)

for i in range(N1):
    zenite_1000GeV[i]=np.arccos(texto1[i][0])
    ratio_1000GeV[i]=texto1[i][1]
func1000GeV=inter.interp1d(zenite_1000GeV,ratio_1000GeV)

theta=np.linspace(3.0*np.pi/180,87*np.pi/180,15)
#print(theta*180/np.pi)
en=np.logspace(0,3,4)
interen=np.zeros(4)
Ei=np.logspace(3,6.69897,1000) #MeV
#print((Ei[0]-np.diff(Ei)[0]/2)/1000)
ratio=np.zeros((len(theta),len(Ei)))
for t in range(len(theta)):
    #print(t)
    interen[0]=func1GeV(theta[t])
    interen[1]=func10GeV(theta[t])
    interen[2]=func100GeV(theta[t])
    interen[3]=func1000GeV(theta[t])
    funcen=inter.interp1d(en,interen)
    for i in range(len(Ei)):
        if (Ei[i]<=1E+6):
            ratio[t][i]=funcen(Ei[i]/1000)
        else:
            ratio[t][i]=2.0*funcen(1E+6/1000)

x=np.linspace(0,(700/np.cos(theta[14]))/700*1855*100,5000)
#print(700*np.diff(x)[0]/1855/100)

E=np.zeros((len(Ei),len(x)))

for i in range(len(Ei)):
    # if (i%500==0):
    #     print(i)
    for j in range(len(x)-1):
        if (j==0):
            E[i][j]=Ei[i]
        
        if (E[i][j]<=1.2): continue
        E[i][j+1]=E[i][j]-funcdEdx(E[i][j])*(x[j+1]-x[j])
        
Acosine=(2*100)**2 #cm^2 COSINE-100 LS Area
Nmu=np.zeros((len(theta),len(Ei)-1))

for i in range(len(Ei)-1):
    # if (i%500==0):
    #     print(i)
    for t in range(len(theta)): 
        I = quad(func, (Ei[i]-np.diff(Ei)[i]/2)/1000, (Ei[i]+np.diff(Ei)[i]/2)/1000)
        sr=2*np.pi*(-np.cos(theta[t]+3.0*np.pi/180)+np.cos(theta[t]-3.0*np.pi/180))
        Nmu[t][i]=ratio[t][i]*I[0]*Acosine*sr #Muons/s

Nmuabove100=np.zeros((len(Ei)-1,len(theta),len(x)))
soma=0.0
for i in range(len(Ei)-1):
    # if (i%500==0):
    #     print(i)
    for t in range(len(theta)):
        j=0
        c=0
        while (E[i][j]>1.2 and np.cos(theta[t])*700*x[j]/1855/100<700):
            if (E[i][j]>1*10**2):
                Nmuabove100[i][t][j]=Nmu[t][i]

            if(np.cos(theta[t])*700*x[j]/1855/100>695 and c==0):
                soma=soma+Nmu[t][i]
                c=c+1
            j=j+1

#print(soma/4*3600*24)

alpha=1/137
me=0.511/1000 #GeV

epsilon=np.logspace(-2,4,30)
totaldecays=np.zeros(len(epsilon))

arq = open('/home/LuisFranca/forDavid/E_z_DPmuonbremss_3gamma.txt', 'r')
texto1 = arq.readlines()
arq.close()

for i in range(len(texto1)):
    texto1[i] = texto1[i].split("\t")
    texto1[i] = [float(x) for x in texto1[i]]

N1=len(texto1)
m_z=np.zeros(N1)
meanE_z=np.zeros(N1)

for i in range(N1):
    m_z[i]=texto1[i][0]
    meanE_z[i]=texto1[i][1]
funcmean=inter.interp1d(m_z,meanE_z)

sig=np.zeros(2000)
arq = open('/home/LuisFranca/forDavid/xsec_DPmuonbremss_m%s_3gamma.txt' % "{:.8f}".format(mDP), 'r')
texto1 = arq.readlines()
arq.close()

for i in range(len(texto1)):
    texto1[i] = texto1[i].split("\t")
    texto1[i] = [float(x) for x in texto1[i]]

N1=len(texto1)
E_mutext=np.zeros(N1)

for i in range(N1):
    E_mutext[i]=texto1[i][0]
    sig[i]=texto1[i][1]
funcsig=inter.interp1d(E_mutext,sig)

f=open('/home/LuisFranca/forDavid/times/ndecays_%s_meanE_3gamma.txt' % mDP,'w')
for k in range(len(epsilon)):
    for i in range(len(Ei)-1):
        if (i%100==0):
            print(i)
        for t in range(len(theta)):
            for j in range(len(x)-1):
                if ((700/np.cos(theta[t])-x[j]/2.65/100)>=0.0):
                    if (E[i][j]<1*10**2): continue
                    xsec=(funcsig(E[i][j]/1000)*(epsilon[k]/1E-4)**2)*10**(-12)*10**(-24)
                    Nbremss=Nmuabove100[i][t][j]*2.65*Na/22*np.diff(x)[0]/2.65*xsec*60*60*24*365.25
                    L=1.973E-16*2*funcmean(mDP)*E[i][j]/1000/(epsilon[k]**2*alpha*mDP**2*(1+2*me**2/mDP**2)*(1-4*me**2/mDP**2)**0.5) # m
                    decays=(np.exp(-((700/np.cos(theta[t])-x[j]/2.65/100)/L))-np.exp(-((702/np.cos(theta[t])-x[j]/2.65/100)/L)))*Nbremss
                    totaldecays[k]=totaldecays[k]+decays
            
    print("{:.3e}".format(epsilon[k]),file=f,end=" ")
    print("{:.3e}".format(totaldecays[k]),file=f,end="\n")
f.close()