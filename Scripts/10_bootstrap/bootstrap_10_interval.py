import itertools, os, os.path, sys, argparse, itertools, shutil
from datetime import datetime
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Bootstrap values for any metric', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('metric', help='metric types, either MCC, ACC, SPE, PRE, F1, or RECALL', choices=['MCC', 'ACC', 'SPE', 'PRE', 'F1', 'RECALL'])
    parser.add_argument('alpha', help='alpha value for significant threshold, either 0.05 or 0.01 is typical',type=float)
    parser.add_argument('reps', help='number of reps used to determine pvalue, either 50 or 100 for 0.05 or 0.01', type=int)
    args = parser.parse_args()


import matplotlib.pyplot as plt
import pandas as pd
import math
import random
import seaborn as sn

import sys
#Make sure this points to the functions folder
#This works assuming you are in analysis
sys.path.append('../Scripts/functions/')

#import all functions
from functions_10_interval import *

#Don't change CARD results
df1 = pd.read_excel('CARD_results.xls')
CARD_dict = {}
rows = []

CARD_dict=mkCARDdict(df1)

#Read in observations
df2 = pd.read_excel('es0c03803_si_002.xls')

phenotype_dict = {}
rows = []
rows=mkrows(df2)
random.shuffle(rows)

phenotype_dict=mkphenofromrow(rows)
df=mkconfusion(CARD_dict,phenotype_dict,pd)


#1 in 10 is 10%
#1 in 50 is 5%
#1 in 100 is 1%
#How many times does the random value
#exceed the observed value
#assigned on command line
alpha=args.alpha
reps=args.reps
#type can be only
#MCC
#ACC
#PRE
#SPE
#RECALL
#F1
metric=args.metric

#for mcc also pass math
#make the obsmat based on chosen metric
if metric=="MCC":
   obsmat=calcmcc(df, math)
elif metric=="ACC":
   obsmat=calcacc(df)
elif metric=="PRE":
   obsmat=calcpre(df)
elif metric=="SPE":
   obsmat=calcspe(df)
elif metric=="F1":
   obsmat=calcf1(df)
elif metric=="RECALL":
   obsmat=calcrecall(df)


maxrmat_dict=mkbootstrap(obsmat, reps, metric, pd, random, math)
# print(maxrmat)
# import pickle
# with open('maxrmat_dict.pickle', 'wb') as handle:
#     pickle.dump(maxrmat_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
# with open('maxrmat_dict.pickle', 'rb') as handle:
#     maxrmat_dict= pickle.load( handle)

# print(maxrmat)
# exit()

current_time = str(datetime.now())
for metric in ["acc","f1","pre","recall","spe","mcc"]:
 maxrmat = maxrmat_dict[metric]
    

 #Figure out whether it's significant or not
 results_mat= [[0 for i in range(11)] for j in range(11)]
 pvalue_mat= [[0 for i in range(11)] for j in range(11)]
 pvalues = []
 for length in range(0,101,10):
  for identity in range(0,101,10):
   #This will give you all pvalues
   #results_mat[length][identity]=maxrmat[length][identity]/reps
   pvalue_mat[length//10][identity//10]=(maxrmat[length//10][identity//10]+1)/(reps+1)
   pvalues.append((maxrmat[length//10][identity//10]+1)/(reps+1))

 # Step 2: Sort the p-values
 pvalues_sorted = sorted(pvalues)
 m = len(pvalues_sorted)

 # Step 3: Apply BH to find the largest i where p_i <= (i/m)*alpha
 threshold_index = -1
 for i in range(m - 1, -1, -1):  # go backwards
     if pvalues_sorted[i] <= alpha * (i + 1) / m:
         threshold_index = i
         break

 # Step 4: Get the BH threshold
 if threshold_index != -1:
     bh_threshold = pvalues_sorted[threshold_index]
 else:
     print("no p-values pass the BH threshold")
     bh_threshold = 0  # no p-values pass the BH threshold

 print(f"{metric} Benjamini-Hochberg threshold: {bh_threshold}")

 # Step 5: Mark results_mat=1 when p-value is NOT significant
 for length in range(0, 11):
     for identity in range(0, 11):
         if pvalue_mat[length][identity] > bh_threshold:
             results_mat[length][identity] = 1



 metrictitle="%s_PVAL Benjamini-Hochberg (alpha %s, blue significant)" % (metric,alpha)
 mkplotRedGreen(results_mat, metrictitle, pd, plt, sn,current_time)
 pfilename="%s/%s_pvalue.txt" % (current_time,metric)
 with open(pfilename, 'w') as f:
  for length in range(0,11,1):
   f.write(f"{pvalue_mat[length]}\n")  
 plt.close()


