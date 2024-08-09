import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.interpolate as inter
import scipy.special as spc
import uproot
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure
import scipy.interpolate as inter

file = uproot.open(f'/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/COSINE Data/Sigma/sigma_DP_SET3_bremss.root')
norm=1E+10
for z in range(0,43):
    for k in range(0,49):
        xsec=[]
        xsec.append(file['hsig%d%d;1' % (z,k)].to_numpy())
        xsec=np.array(xsec, dtype="object")
        xsec[0][0]=xsec[0][0]/np.max(xsec[0][0])
        binxsecmed=np.zeros(len(xsec[0][1])-1)
        newxsecmed=np.zeros(51)
        xsecnewxsec=np.zeros(50)
        mean=0.0
        counter=0
        for i in range(len(xsec[0][1])-1):
            if (i!=len(xsec[0][1])-1):
                binxsecmed[i]=(xsec[0][1][i]+xsec[0][1][i+1])/2
                mean=mean+xsec[0][0][i]
            if (i%20==0 and i!=0):
                newxsecmed[counter+1]=binxsecmed[i]
                xsecnewxsec[counter]=mean/20
                mean=0.0
                counter=counter+1
            newxsecmed[50]=binxsecmed[999]

        integr=0
        for i in range(len(binxsecmed)-1):
            integr=integr+(xsec[0][0][i]+xsec[0][0][i+1])*(binxsecmed[i+1]-binxsecmed[i])/2
        cdf=np.zeros(len(binxsecmed))
        totalint=integr
        lowen=binxsecmed[0]
        highen=lowen+np.diff(binxsecmed)[0]
        for i in range(len(binxsecmed)-1):
            integr=(xsec[0][0][i]+xsec[0][0][i+1])*(binxsecmed[i+1]-binxsecmed[i])/2
            if (i!=0):
                cdf[i+1]=cdf[i]+integr/totalint
            else:
                cdf[i]=integr/totalint
        cdffunc=inter.interp1d(cdf,binxsecmed, fill_value="extrapolate")
        upper_data_SET3=cdffunc(0.9)*norm
        print(upper_data_SET3)

        f1=plt.figure()
        plt.title("Posterior for DP event excess")
        plt.stairs(xsecnewxsec*1.0/np.max(xsecnewxsec),newxsecmed*norm,linewidth=2.0,color='black',label='Scaled posterior PDF')
        plt.plot(binxsecmed*norm,cdf,linewidth=1.5,label='CDF',color='green')
        plt.vlines(upper_data_SET3,0,1,linestyle='dashed',color='red',linewidth=1.3,label='90% CL')
        #plt.xlim(0,3)
        plt.ylim(0,1.0)
        #plt.xlabel('$\sigma_e$ ($10^{-31}\;cm^2$)',fontsize=15)
        plt.ylabel('Probability',fontsize=12)
        plt.xlabel('Signal strength', fontsize=12)
        plt.legend(fontsize=11,loc='center right')
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        #plt.show()
        plt.savefig("/home/davidc/Documents/Master's Analisys/Parameter space/Codes/COSINE Data/Posteriors/Posterior_%dm_%dk.png" %(z, k))
        plt.close()