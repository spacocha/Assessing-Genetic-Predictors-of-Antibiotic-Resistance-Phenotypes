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

#make sure your 1 path (don't use 0) points to functions folder
sys.path.append('../Assessing-Genetic-Predictors-of-Antibiotic-Resistance-Phenotypes/Scripts/functions/')

#import all functions
from functions import *

sys.path.append('../Assessing-Genetic-Predictors-of-Antibiotic-Resistance-Phenotypes/Scripts/08_confusion_matrix/')
from confusion_matrix import df

#df comes from confusion script
matrix=calcmcc(df,math)
type='MCC'
mkplot(matrix, type, pd, plt, sn)
plt.close()

#run accuracy
matrix=calcacc(df)
type='ACC'
mkplot(matrix, type, pd, plt, sn)
plt.close()

matrix=calcpre(df)
type='PRE'
mkplot(matrix, type, pd, plt, sn)
plt.close()

type='SPE'
matrix=calcspe(df)
mkplot(matrix, type, pd, plt, sn)
plt.close()

type='RECALL'
matrix=calcrecall(df)
mkplot(matrix, type, pd, plt, sn)
plt.close()

type='F1'
matrix=calcf1(df)
mkplot(matrix, type, pd, plt, sn)
plt.close()


