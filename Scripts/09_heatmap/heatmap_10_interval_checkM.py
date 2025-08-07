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
from functions_10_interval_checkM import *

sys.path.append('../Scripts/08_confusion_matrix/')
from confusion_matrix_10_interval_checkM import df


matrix=calcpre(df)
type='PRE after checkM'
mkplot(matrix, type, pd, plt, sn)
plt.close()
