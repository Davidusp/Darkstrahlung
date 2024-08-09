from cProfile import label
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import scipy.interpolate as inter
import pandas as pd
import sys
import os

def main():

    mDP=float(sys.argv[1])
    epsilon=float(sys.argv[2])
    # mDP=0.0011514
    # epsilon=0.00000000010
    
    Ebins = np.linspace(9, 1485, 739)
    Binedges = np.linspace(8, 1486, 740) #binsize = 2 NPE

    arq2 = open('/data/COSINE/WORK/dfreitas/single_hit/CorrectedSpec100keV/CorrectedSpectrum_m=%.8f_k=%.11f.txt' %(mDP, epsilon ),'r')
    #arq2 = open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/COSINE Data/CorrectedSpectra/CorrectedSpectrum_m=%.8f_k=%.11f.txt' %(mDP, epsilon ),'r')

    NPEC2 = []
    NPEC3 = []
    NPEC4 = []
    NPEC6 = []
    NPEC7 = []

    SpecC2 = []
    SpecC3 = []
    SpecC4 = []
    SpecC6 = []
    SpecC7 = []

    SpecnewC2 = []
    SpecnewC3 = []
    SpecnewC4 = []
    SpecnewC6 = []
    SpecnewC7 = []

    for linha in arq2:
        line = linha.split("	")
        NPEC2.append(float(line[0]))
        SpecC2.append(float(line[1]))
        NPEC3.append(float(line[2]))
        SpecC3.append(float(line[3]))
        NPEC4.append(float(line[4]))
        SpecC4.append(float(line[5]))
        NPEC6.append(float(line[6]))
        SpecC6.append(float(line[7]))
        NPEC7.append(float(line[8]))
        SpecC7.append(float(line[9]))

    arq2.close()

    sumC2 = 0
    sumC3 = 0
    sumC4 = 0
    sumC6 = 0
    sumC7 = 0

    countC2 = 1
    countC3 = 1
    countC4 = 1
    countC6 = 1
    countC7 = 1

    posC2 = 0
    posC3 = 0
    posC4 = 0
    posC6 = 0
    posC7 = 0
    
    for i in range(0,len(NPEC2)-1):
        if NPEC2[i] >= Binedges[0] and NPEC2[i] <= Binedges[len(Binedges)-1]:
            if NPEC2[i] >= Binedges[countC2-1] and NPEC2[i] < Binedges[countC2]:
                sumC2 += SpecC2[i]

            else:
                SpecnewC2.append(sumC2/(i+1-posC2)) #Computing the mean values in each energy bin
                sumC2 = SpecC2[i]
                posC2 = i+1
                countC2 += 1
        
        elif i>0:
            if NPEC2[i-1] >= Binedges[0] and NPEC2[i-1] <= Binedges[len(Binedges)-1]:
                SpecnewC2.append(sumC2/(i+1-posC2)) #Computing the mean values in each energy bin

        if NPEC3[i] >= Binedges[0] and NPEC3[i] <= Binedges[len(Binedges)-1]:
            if NPEC3[i] >= Binedges[countC3-1] and NPEC3[i] < Binedges[countC3]:
                sumC3 += SpecC3[i]
            else:
                SpecnewC3.append(sumC3/(i+1-posC3)) #Computing the mean values in each energy bin
                sumC3 = SpecC3[i]
                posC3 = i+1
                countC3 += 1
        elif i>0:
            if NPEC3[i-1] >= Binedges[0] and NPEC3[i-1] <= Binedges[len(Binedges)-1]:
                SpecnewC3.append(sumC3/(i+1-posC3)) #Computing the mean values in each energy bin

        if NPEC4[i] >= Binedges[0] and NPEC4[i] <= Binedges[len(Binedges)-1]:
            if NPEC4[i] >= Binedges[countC4-1] and NPEC4[i] < Binedges[countC4]:
                sumC4 += SpecC4[i]
            else:
                SpecnewC4.append(sumC4/(i+1-posC4)) #Computing the mean values in each energy bin
                sumC4 = SpecC4[i]
                posC4 = i+1
                countC4 += 1
        elif i>0:
            if NPEC4[i-1] >= Binedges[0] and NPEC4[i-1] <= Binedges[len(Binedges)-1]:
                SpecnewC4.append(sumC4/(i+1-posC4)) #Computing the mean values in each energy bin

        if NPEC6[i] >= Binedges[0] and NPEC6[i] <= Binedges[len(Binedges)-1]:
            if NPEC6[i] >= Binedges[countC6-1] and NPEC6[i] < Binedges[countC6]:
                sumC6 += SpecC6[i]
            else:
                SpecnewC6.append(sumC6/(i+1-posC6)) #Computing the mean values in each energy bin
                sumC6 = SpecC6[i]
                posC6 = i+1
                countC6 += 1
        elif i>0:
            if NPEC6[i-1] >= Binedges[0] and NPEC6[i-1] <= Binedges[len(Binedges)-1]:
                SpecnewC6.append(sumC6/(i+1-posC6)) #Computing the mean values in each energy bin
        
        if NPEC7[i] >= Binedges[0] and NPEC7[i] <= Binedges[len(Binedges)-1]:
            if NPEC7[i] >= Binedges[countC7-1] and NPEC7[i] < Binedges[countC7]:
                sumC7 += SpecC7[i]
            else:
                SpecnewC7.append(sumC7/(i+1-posC7)) #Computing the mean values in each energy bin
                sumC7 = SpecC7[i]
                posC7 = i+1
                countC7 += 1
        elif i>0:
            if NPEC7[i-1] >= Binedges[0] and NPEC7[i-1] <= Binedges[len(Binedges)-1]:
                SpecnewC7.append(sumC7/(i+1-posC7)) #Computing the mean values in each energy bin

    # print(len(SpecnewC2),len(SpecnewC3),len(SpecnewC4),len(SpecnewC6),len(SpecnewC7))
    # print(countC2,countC3,countC4,countC6,countC7)
    # form_mDP = f"{mDP:.3e}"
    # form_kappa = f"{epsilon:.3e}"
    # file_name = f"Specorr100keV_m={form_mDP}_k={form_kappa}.txt"
    # file_path = os.path.join("/data/COSINE/WORK/dfreitas/single_hit/FinalSpec100keV", file_name)
    #file_path = os.path.join("/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/COSINE Data", file_name)

    with open('/data/COSINE/WORK/dfreitas/single_hit/FinalSpec100keV/Specorr100keV_m=%d_k=%d.txt' %(round(mDP*1e07),round(epsilon*1e11)), 'w') as archive:
        for i in range(0, len(Ebins)):
            archive.write('{:.1f}    {:.16f}    {:.16f}    {:.16f}    {:.16f}    {:.16f}\n'.format(Ebins[i], SpecnewC2[i], SpecnewC3[i], SpecnewC4[i], SpecnewC6[i], SpecnewC7[i]))

    # with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/COSINE Data/FinalSpectra/Specorr100keV_m=%d_k=%d.txt' %(round(mDP*1e07),round(epsilon*1e11)), 'w') as archive:
    #     for i in range(0, len(Ebins)):
    #         archive.write('{:.1f}    {:.16f}    {:.16f}    {:.16f}    {:.16f}    {:.16f}\n'.format(Ebins[i], SpecnewC2[i], SpecnewC3[i], SpecnewC4[i], SpecnewC6[i], SpecnewC7[i]))
            
main()