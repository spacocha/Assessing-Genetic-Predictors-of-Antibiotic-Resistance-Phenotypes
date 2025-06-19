import itertools, os, os.path, sys, argparse, itertools, shutil

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
from functions import *

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


maxrmat=mkbootstrap(obsmat, reps, metric, pd, random, math)

#Figure out whether it's significant or not
results_mat= [[0 for i in range(101)] for j in range(101)]
pvalue_mat= [[0 for i in range(101)] for j in range(101)]
for length in range(0,101,1):
 for identity in range(0,101,1):
  #This will give you all pvalues
  #Pvalue is more accurately rep as (r+1)/(n+1)
  #Am. J. Hum. Genet. 71:439–441, 2002
  pvalue_mat[length][identity]=(maxrmat[length][identity]+1)/(reps+1)
  if pvalue_mat[length][identity] > alpha:
   results_mat[length][identity]=1

pfilename="%s_pvalue.txt" % metric
with open(pfilename, 'w') as f:
 for length in range(0,101,1):
  f.write(f"{pvalue_mat[length]}\n")  

metrictitle="%s_PVAL" % metric
mkplot(results_mat, metrictitle, pd, plt, sn)
plt.close()


