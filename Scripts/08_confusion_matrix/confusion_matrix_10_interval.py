#might have to do this eventually
#Right now, put functions.py where you execute this
#
import sys

#make sure this points to functions folder
#This will work assuming you are working in analysis
sys.path.append('../Scripts/functions/')

#import all functions
from functions_10_interval import *

import pandas as pd
import math
import random

#not sure what this small test is
#xls = pd.ExcelFile("small data set test.xlsx")
#changed to import excel files directly
df1 = pd.read_excel('CARD_results.xls')
#remove sheet 2 from Gray et al SI table 2
df2 = pd.read_excel('es0c03803_si_002.xls')

#old phenotype data by class
#df2 = pd.read_excel('phenotype data.xls')

phenotype_dict = {}
rows = []
rows=mkrows(df2)
random.shuffle(rows)

phenotype_dict=mkphenofromrow(rows)

#Is this used?
#resistance_dict = {}

CARD_dict = {}
rows = []
CARD_dict=mkCARDdict(df1)

#Is this used? If so, where?
import seaborn as sn

df=mkconfusion(CARD_dict,phenotype_dict,pd)
