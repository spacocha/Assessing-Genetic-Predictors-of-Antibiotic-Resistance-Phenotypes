import math
import matplotlib.pyplot as plt
import seaborn as sn
import pandas as pd
import random

# continue from 08_confusion_matrix
#import functions from functions.py
#for now, just copy/past
#I need to figure out how to keep df from the confusion script

import sys

#make sure this points to functions folder
#This works assuming you are in analysis
sys.path.append('../Scripts/functions/')

#import all functions
from functions_10_interval_checkM_confusion import *

#Read in raw data
df1 = pd.read_csv('CARD_results.csv')
df1 = df1.fillna('')
df2 = pd.read_csv('es0c03803_si_002.csv')
df3 = pd.read_csv('checkm2_results.csv')

#make confusion matrix
df=confusion_matrix_10_interval_checkM_def(df1, df2, df3, pd, math, random)

#Calculate the metrics
type="PRE_after_checkM"
with open('%s_df_output.txt' %type, '+w') as f:
 for items in df:
  f.write('%s\n' %items)

matrix=calcpre(df)

with open('%s_matrix_output.txt' %type, '+w') as f:
 for items in matrix:
  f.write('%s\n' %items)

mkplot(matrix, type, pd, plt, sn)
plt.close()

type="mcc_after_checkM"
matrix=calcmcc(df,math)

with open('%s_matrix_output.txt' %type, '+w') as f:
 for items in matrix:
  f.write('%s\n' %items)

mkplot(matrix, type, pd, plt, sn)
plt.close()

type="acc_after_checkM"
matrix=calcacc(df)

with open('%s_matrix_output.txt' %type, '+w') as f:
 for items in matrix:
  f.write('%s\n' %items)

mkplot(matrix, type, pd, plt, sn)
plt.close()

type="RECALL_after_checkM"
matrix=calcrecall(df)

with open('%s_matrix_output.txt' %type, '+w') as f:
 for items in matrix:
  f.write('%s\n' %items)

mkplot(matrix, type, pd, plt, sn)
plt.close()

type="SPE_after_checkM"
matrix=calcspe(df)

with open('%s_matrix_output.txt' %type, '+w') as f:
 for items in matrix:
  f.write('%s\n' %items)

mkplot(matrix, type, pd, plt, sn)
plt.close()

type="F1_after_checkM"
matrix=calcf1(df)

with open('%s_matrix_output.txt' %type, '+w') as f:
 for items in matrix:
  f.write('%s\n' %items)

mkplot(matrix, type, pd, plt, sn)
plt.close()
