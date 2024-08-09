import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.interpolate as inter
import sys


mDP = np.array([1,2,3,4,5,7,10,13,17,22,29,39,52,69,91,121,160,212,281,373,494,655,869,1151,1526,11514,2024,
                2683,3556,4715,6251,8286,10985,14563,19307,25595,33932,44984,59636,79060,104811,138950,184207])

kappas = np.array([10,15,24,37,56,87,133,205,316,487,750,1155,1778,2738,4217,6494,10000,15399,23714,36517,56234,
                   86596,133352,205353,316228,486968,749894,1154782,1778279,2738420,4216965,6493816,10000000,15399265,
                   23713737,36517413,56234133,86596432,133352143,205352503,316227766,486967525,749894209,1154781985,
                   1778279410,2738419634,4216965034,6493816316,10000000000])

for z in range(0,len(mDP)):
    for k in range(0,len(kappas)):
        with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/COSINE Data/FinalSpectra/Specorr100keV_m=%d_k=%d.txt' %(mDP[z],kappas[k]), 'r') as archive:
            print("cool")