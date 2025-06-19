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
from functions import *

sys.path.append('../Scripts/08_confusion_matrix/')
from confusion_matrix import df

#df comes from confusion script
matrix=calcmcc(df,math)
metric='MCC'
mkplot(matrix, metric, pd, plt, sn)
plt.close()

#run accuracy
matrix=calcacc(df)
metric='ACC'
mkplot(matrix, metric, pd, plt, sn)
plt.close()

matrix=calcpre(df)
metric='PRE'
mkplot(matrix, metric, pd, plt, sn)
plt.close()

metric='SPE'
matrix=calcspe(df)
mkplot(matrix, metric, pd, plt, sn)
plt.close()

metric='RECALL'
matrix=calcrecall(df)
mkplot(matrix, metric, pd, plt, sn)
plt.close()

metric='F1'
matrix=calcf1(df)
mkplot(matrix, metric, pd, plt, sn)
plt.close()


