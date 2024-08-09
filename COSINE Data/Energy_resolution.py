import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.interpolate as inter
import scipy.special as spc
import sys
import uproot

# Mass and day values given by user. It makes this script parallelizable.
mdm=float(sys.argv[1])   
kappa = float(sys.argv[2])
# mdm = 0.00202400   
# kappa = 0.00056234133
me=0.511
 
energy=np.zeros(100) #Energies (100)
spectramax_Na_v=np.zeros(len(energy)) #Detection events
                
#arq = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/CompleteSpec/SpecAllEarth/100keV/NDir100keV_int_m=%.8f_k=%.11f.txt' %(mdm,kappa), 'r')
arq = open('/data/COSINE/WORK/dfreitas/single_hit/Spec100keVFR/NDir100keV_int_m=%.8f_k=%.11f.txt' %(mdm,kappa), 'r')
arq2 = open('/data/COSINE/WORK/dfreitas/single_hit/SpecAbv100keV/NDir100keV_int_m=%.8f_k=%.11f.txt' %(mdm,kappa),'r')
texto1 = arq.readlines()
texto2 = arq2.readlines()
arq.close()
arq2.close()

for i in range(len(texto1)):
    texto1[i] = texto1[i].split("  ")
    texto2[i] = texto2[i].split("  ")
    texto1[i] = [float(x) for x in texto1[i]]
    texto2[i] = [float(y) for y in texto2[i]]


N1=100
for i in range(0,N1):
    energy[i]=texto1[i][0]*1000 #keV
    #if mdm[i]<2*me:
    spectramax_Na_v[i]=(texto1[i][1]+texto1[i][2]+texto2[i][1]+texto2[i][2]+texto2[i][3])/106/365.25 #events/kg/day
    #print(spectramax_Na_v[i], "spectrum")
    # else:
    #     spectramax_Na_v[i]=texto1[i][1]/106/365.25

#file = uproot.open(f'/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/COSINE Data/NPE_SET3_nPR.root')
file = uproot.open(f'/data/COSINE/WORK/dfreitas/single_hit/NPE_SET3_nPR.root')
energyc2=file['npr_graph_2'].member('fX')
ratioc2=file['npr_graph_2'].member('fY')
entoratio_c2func=inter.interp1d(energyc2,ratioc2)
energyc3=file['npr_graph_3'].member('fX')
ratioc3=file['npr_graph_3'].member('fY')
entoratio_c3func=inter.interp1d(energyc3,ratioc3)
energyc4=file['npr_graph_4'].member('fX')
ratioc4=file['npr_graph_4'].member('fY')
entoratio_c4func=inter.interp1d(energyc4,ratioc4)
energyc6=file['npr_graph_6'].member('fX')
ratioc6=file['npr_graph_6'].member('fY')
entoratio_c6func=inter.interp1d(energyc6,ratioc6)
energyc7=file['npr_graph_7'].member('fX')
ratioc7=file['npr_graph_7'].member('fY')
entoratio_c7func=inter.interp1d(energyc7,ratioc7)

qc550keVc2=747.3*620.697571
qc550keVc3=747.8*639.862549
qc550keVc4=750*587.43811
qc550keVc6=745.3*565.981445
qc550keVc7=748.2*590.684814

qc5c2=np.zeros(100)
qc5c3=np.zeros(100)
qc5c4=np.zeros(100)
qc5c6=np.zeros(100)
qc5c7=np.zeros(100)

for i in range(len(energy)):
    if energy[i] > 20000: continue
    if (energy[i]>=0.05 and energy[i]<200):
        qc5c2[i]=entoratio_c2func(energy[i])*qc550keVc2/50.0*energy[i]
        qc5c3[i]=entoratio_c3func(energy[i])*qc550keVc3/50.0*energy[i]
        qc5c4[i]=entoratio_c4func(energy[i])*qc550keVc4/50.0*energy[i]
        qc5c6[i]=entoratio_c6func(energy[i])*qc550keVc6/50.0*energy[i]
        qc5c7[i]=entoratio_c7func(energy[i])*qc550keVc7/50.0*energy[i]
    elif (energy[i]>=200 and energy[i]<20000):
        qc5c2[i]=1.0*qc550keVc2/50.0*energy[i]
        qc5c3[i]=1.0*qc550keVc3/50.0*energy[i]
        qc5c4[i]=1.0*qc550keVc4/50.0*energy[i]
        qc5c6[i]=1.0*qc550keVc6/50.0*energy[i]
        qc5c7[i]=1.0*qc550keVc7/50.0*energy[i]
        
f2=np.zeros(len(qc5c2))
f3=np.zeros(len(qc5c3))
f4=np.zeros(len(qc5c4))
f6=np.zeros(len(qc5c6))
f7=np.zeros(len(qc5c7))

#p0=[0.001789,0.001885,0.000724,0.001725,0.001258]
#p1=[918.3219,785.8450,968.3337,694.9439,825.2153]
p0=[0.00068729,0.000584246,0.00057917,0.00078596,0.00043378]
p1=[888.12252,844.5472,901.56205,771.9161,917.1382]
for i in range(len(qc5c2)):
    
    if (qc5c2[i]!=0):
        f2[i]=(p0[0]+p1[0]/qc5c2[i])**0.5*qc5c2[i]
    if (qc5c3[i]!=0):
        f3[i]=(p0[1]+p1[1]/qc5c3[i])**0.5*qc5c3[i]
    if (qc5c4[i]!=0):
        f4[i]=(p0[2]+p1[2]/qc5c4[i])**0.5*qc5c4[i]
    if (qc5c6[i]!=0):
        f6[i]=(p0[3]+p1[3]/qc5c6[i])**0.5*qc5c6[i]
    if (qc5c7[i]!=0):
        f7[i]=(p0[4]+p1[4]/qc5c7[i])**0.5*qc5c7[i]
    
qapp_gauss=np.linspace(0,2e06,20000)
G2=np.zeros((len(qc5c2),len(qapp_gauss)))
G3=np.zeros((len(qc5c3),len(qapp_gauss)))
G4=np.zeros((len(qc5c4),len(qapp_gauss)))
G6=np.zeros((len(qc5c6),len(qapp_gauss)))
G7=np.zeros((len(qc5c7),len(qapp_gauss)))

for j in range(len(qc5c2)):
    print(qc5c2[j],spectramax_Na_v[j])
    #if (qc5c2[j]<620.697571): continue
    if qc5c2[j]!=0:
        s2=f2[j]**2/qc5c2[j]
        integrG2=0.0
        for i in range(len(qapp_gauss)):
            if (np.abs(qc5c2[j]-qapp_gauss[i])>20000000): continue
            G2[j][i]=np.exp(-0.5*((qapp_gauss[i]-qc5c2[j])/f2[j])**2)
            if (i!=0):
                integrG2=integrG2+(G2[j][i]+G2[j][i-1])*(qapp_gauss[i]-qapp_gauss[i-1])/2
        G2[j]=G2[j]/integrG2

        #if j%(len(qc5c2)-1)==0:
        #if j==9:
        #    print(G2[j], qc5c2[j])

        #    f1=plt.figure()
        #    plt.plot(qapp_gauss, G2[j])
        #    plt.xlabel("qapp_gauss", fontsize=11)
        #    plt.ylabel("$G2[0]$", fontsize=12)
        #    plt.show()
        #    plt.close()
    
#exit()
for j in range(len(qc5c3)):
    #if (qc5c3[j]<639.862549): continue
    if qc5c3[j]!=0:
        s3=f3[j]**2/qc5c3[j]
        integrG3=0.0
        for i in range(len(qapp_gauss)):
            G3[j][i]=np.exp(-0.5*((qapp_gauss[i]-qc5c3[j])/f3[j])**2)
            if (i!=0):
                integrG3=integrG3+(G3[j][i]+G3[j][i-1])*(qapp_gauss[i]-qapp_gauss[i-1])/2
        G3[j]=G3[j]/integrG3

for j in range(len(qc5c4)):
    #if (qc5c4[j]<587.43811): continue
    if qc5c4[j]!=0:
        s4=f4[j]**2/qc5c4[j]
        integrG4=0.0
        for i in range(len(qapp_gauss)):
            G4[j][i]=np.exp(-0.5*((qapp_gauss[i]-qc5c4[j])/f4[j])**2)
            if (i!=0):
                integrG4=integrG4+(G4[j][i]+G4[j][i-1])*(qapp_gauss[i]-qapp_gauss[i-1])/2
        G4[j]=G4[j]/integrG4

for j in range(len(qc5c6)):
    #if (qc5c6[j]<565.981445): continue
    if qc5c6[j]!=0:
        s6=f6[j]**2/qc5c6[j]
        integrG6=0.0
        integrS6=0.0
        for i in range(len(qapp_gauss)):
            G6[j][i]=np.exp(-0.5*((qapp_gauss[i]-qc5c6[j])/f6[j])**2)
            if (i!=0):
                integrG6=integrG6+(G6[j][i]+G6[j][i-1])*(qapp_gauss[i]-qapp_gauss[i-1])/2
        G6[j]=G6[j]/integrG6

for j in range(len(qc5c7)):
    #if (qc5c7[j]<590.684814): continue
    if qc5c7[j]!=0:
        s7=f7[j]**2/qc5c7[j]
        integrG7=0.0
        for i in range(len(qapp_gauss)):
            G7[j][i]=np.exp(-0.5*((qapp_gauss[i]-qc5c7[j])/f7[j])**2)
            if (i!=0):
                integrG7=integrG7+(G7[j][i]+G7[j][i-1])*(qapp_gauss[i]-qapp_gauss[i-1])/2
        G7[j]=G7[j]/integrG7
    
raterescorr_c2_gauss=np.zeros(len(qapp_gauss))
raterescorr_c3_gauss=np.zeros(len(qapp_gauss))
raterescorr_c4_gauss=np.zeros(len(qapp_gauss))
raterescorr_c6_gauss=np.zeros(len(qapp_gauss))
raterescorr_c7_gauss=np.zeros(len(qapp_gauss))
integrc2_gauss=np.zeros(len(qc5c2))
integrc3_gauss=np.zeros(len(qc5c3))
integrc4_gauss=np.zeros(len(qc5c4))
integrc6_gauss=np.zeros(len(qc5c6))
integrc7_gauss=np.zeros(len(qc5c7))
for j in range(len(qapp_gauss)):
    for i in range(len(qc5c2)):
        integrc2_gauss[i]=spectramax_Na_v[i]*G2[i][j]
        integrc3_gauss[i]=spectramax_Na_v[i]*G3[i][j]
        integrc4_gauss[i]=spectramax_Na_v[i]*G4[i][j]
        integrc6_gauss[i]=spectramax_Na_v[i]*G6[i][j]
        integrc7_gauss[i]=spectramax_Na_v[i]*G7[i][j]
        if (i!=0 and qc5c2[i]!=0):
            raterescorr_c2_gauss[j]=raterescorr_c2_gauss[j]+(integrc2_gauss[i]+integrc2_gauss[i-1])*(qc5c2[i]-qc5c2[i-1])/2
        if (i!=0 and qc5c3[i]!=0):
            raterescorr_c3_gauss[j]=raterescorr_c3_gauss[j]+(integrc3_gauss[i]+integrc3_gauss[i-1])*(qc5c3[i]-qc5c3[i-1])/2
        if (i!=0 and qc5c4[i]!=0):
            raterescorr_c4_gauss[j]=raterescorr_c4_gauss[j]+(integrc4_gauss[i]+integrc4_gauss[i-1])*(qc5c4[i]-qc5c4[i-1])/2
        if (i!=0 and qc5c6[i]!=0):
            raterescorr_c6_gauss[j]=raterescorr_c6_gauss[j]+(integrc6_gauss[i]+integrc6_gauss[i-1])*(qc5c6[i]-qc5c6[i-1])/2
        if (i!=0 and qc5c7[i]!=0):
            raterescorr_c7_gauss[j]=raterescorr_c7_gauss[j]+(integrc7_gauss[i]+integrc7_gauss[i-1])*(qc5c7[i]-qc5c7[i-1])/2

# plt.plot(qapp_gauss, raterescorr_c2_gauss)
# plt.show()

f2=open('/data/COSINE/WORK/dfreitas/single_hit/CorrectedSpec100keV/CorrectedSpectrum_m=%.8f_k=%.11f.txt' % (mdm, kappa),'w')
for j in range(len(qapp_gauss)):
    print("{:.5e}".format(qapp_gauss[j]/620.697571),file=f2,end="\t")
    print("{:.5e}".format(raterescorr_c2_gauss[j]),file=f2,end="\t")
    print("{:.5e}".format(qapp_gauss[j]/639.862549),file=f2,end="\t")
    print("{:.5e}".format(raterescorr_c3_gauss[j]),file=f2,end="\t")
    print("{:.5e}".format(qapp_gauss[j]/587.43811),file=f2,end="\t")
    print("{:.5e}".format(raterescorr_c4_gauss[j]),file=f2,end="\t")
    print("{:.5e}".format(qapp_gauss[j]/565.981445),file=f2,end="\t")
    print("{:.5e}".format(raterescorr_c6_gauss[j]),file=f2,end="\t")
    print("{:.5e}".format(qapp_gauss[j]/590.684814),file=f2,end="\t")
    print("{:.5e}".format(raterescorr_c7_gauss[j]),file=f2,end="\n")
f2.close()
