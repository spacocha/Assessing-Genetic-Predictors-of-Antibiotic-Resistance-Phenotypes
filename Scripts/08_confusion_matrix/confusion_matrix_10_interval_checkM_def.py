#might have to do this eventually
#Right now, put functions.py where you execute this
#
import sys

#make sure this points to functions folder
#This will work assuming you are working in analysis
sys.path.append('../Scripts/functions/')

#import all functions
from functions_10_interval_checkM_confusion import *

import pandas as pd
import math
import random

df1 = pd.read_csv('CARD_results.csv')
df1 = df1.fillna('')
df2 = pd.read_csv('es0c03803_si_002.csv')
df3 = pd.read_csv('checkm2_results.csv')

df=confusion_matrix_10_interval_checkM_def(df1, df2, df3, pd, math, random)

