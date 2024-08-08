import math
from math import isnan, nan
from os import O_SYNC, kill
from tracemalloc import stop
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy.interpolate as inter
import sys
import numba

def main():
    with open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/Direct Interaction/NDir_int_m=0.05963600_k=0.01778279410.txt', 'r') as archive:
            arq = archive.read()

    line= arq.split("\n")
    line = line[:-1]
    
    for row in range(0, len(line)):
        line[row] = line[row].split("  ")
        line[row] = [float(x) for x in line[row]]

    events = 0

    for l in range(0, len(line)):
        events += line[l][1] + line[l][2] + line[l][3]

    print(events)
        
main()