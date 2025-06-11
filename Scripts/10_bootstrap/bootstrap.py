import matplotlib.pyplot as plt
import pandas as pd
import math
import random

import sys
sys.path.append('../Assessing-Genetic-Predictors-of-Antibiotic-Resistance-Phenotypes/Scripts/functions/')

#import all functions
from functions import *

#Don't change CARD results
df1 = pd.read_excel('CARD results.xls')
CARD_dict = {}
rows = []

CARD_dict=mkCARDdict(df1)


#1 in 10 is 10%
#1 in 50 is 5%
#1 in 100 is 1%
#This doesn't account for multiple
#hypothesis testing
maxrmat=[]
for i in range(0,50):
 #Create 50 random matrices in the program
 #Read in real phenotypes
 #remove sheet 2 from Gray et al SI table 2
 df2 = pd.read_excel('es0c03803_si_002.xls')
 #shuffle each of the important antibiotics
 random.shuffle(df2.amp_res)
 random.shuffle(df2.cip_res)
 random.shuffle(df2.tet_res)
 random.shuffle(df2.c_res)
 random.shuffle(df2.gm_res)
 random.shuffle(df2.azm_res)
 random.shuffle(df2.cl_res)
 #clear variables
 phenotype_dict = {}
 rows = []
 #remake the phenotype_dict
 rows=mkrows(df2)
 random.shuffle(rows)
 phenotype_dict=mkphenofromrow(rows)
 #Is this used?
 import seaborn as sn
 #make confusion matrix
 df=mkconfusion(CARD_dict,phenotype_dict,pd)
 #calculate acc
 matrix=calcpre(df)
 #max value of random matrices compared to real values
 if i==0:
  #First, so set maxrmat to matrix
  maxrmat=matrix
 else:
  for length in range(0,121,1):
   for identity in range(0,101,1):
    if matrix[length][identity]>maxrmat[length][identity]:
     maxrmat[length][identity]=matrix[length][identity]


#make the true observations matrix
df2 = pd.read_excel('es0c03803_si_002.xls')
phenotype_dict = {}
rows = []
rows=mkrows(df2)
phenotype_dict=mkphenofromrow(rows)
df=mkconfusion(CARD_dict,phenotype_dict,pd)

#for mcc also pass math
obsmat=calcpre(df)

import seaborn as sn

result_matrix = [[0 for i in range(101)] for j in range(121)]

#pvalue 0 if < 0.01
#i.e. value this high not see by chance
#pvalue 1 if > 0.01
#i.e. value this high see by chance 1% of time
for length in range(0,121,1):
 for identity in range(0,101,1):
  if maxrmat[length][identity] < obsmat[length][identity]:
   result_matrix[length][identity] =0
  else:
   result_matrix[length][identity] = 1

type='PRE_PVAL'
mkplot(result_matrix, type, pd, plt, sn)
plt.close()

#Add the other metrics, or save as separate files
