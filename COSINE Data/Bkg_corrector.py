from scipy.optimize import curve_fit
import uproot
import numpy as np
import glob

infiles=[file+':t' for file in glob.glob(f'/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/COSINE Data/FitResult_C6.root')]
activity = uproot.concatenate(infiles, filter_name=['act'],library='numpy')
activity=np.array(activity['act'])

file = uproot.open(f'/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/COSINE Data/FitResult_C6.root')
f1=open('/home/davidc/Documents/Master\'s Analisys/Parameter space/Codes/COSINE Data/FinalBkg/SET3_bkg_C6_8xpe_poisson_nonpr.txt','w')
xpe=np.linspace(9, 1485, 739)
allbkg=np.zeros((57,739))
lenbin=739
surfpb210=np.zeros(lenbin)
tes=np.zeros(lenbin)
#others=np.zeros(lenbin)
for component in range(57):
    bkg=file['hRslt%s_1_1;1' % (component+1)].to_numpy()
    for j in range(lenbin):
        allbkg[component][j]=bkg[0][j]

#First 7 components
indiv_comp = [3,4,10,11,12,13,25]

#other components
others = [1,2,5,14,15,16,17,18,19,20,21,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57]

for j in range(lenbin):
    sum = 0
    for component in range(57):

        #if (rate_poisson[component][j]==0.0): rate_poisson[component][j]=1E-20
        if (component+1 in indiv_comp):
            print("{:.5e}".format(allbkg[component][j]),file=f1,end=" ")

        if ((component+1)==6):  
            print("{:.5e}".format(allbkg[5][j]+allbkg[6][j]+allbkg[7][j]+allbkg[8][j]),file=f1,end=" ")
           
        if ((component+1)==22):
            print("{:.5e}".format(allbkg[21][j]+allbkg[22][j]+allbkg[23][j]),file=f1,end=" ")
        
        if (component in others):
            sum += allbkg[component-1][j]
        
        # if (component!=(56)):
        #     print("{:.5e}".format(sum),file=f1,end=" ")
        if (component==(56)):
            print("{:.5e}".format(sum),file=f1,end="\n")
        
f1.close()