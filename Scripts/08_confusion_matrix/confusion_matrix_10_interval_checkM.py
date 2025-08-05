#might have to do this eventually
#Right now, put functions.py where you execute this
#
import sys

#make sure this points to functions folder
#This will work assuming you are working in analysis
sys.path.append('../Scripts/functions/')

#import all functions
from functions_10_interval_checkM import *

import pandas as pd
import math
import random

#not sure what this small test is
#xls = pd.ExcelFile("small data set test.xlsx")
#changed to import excel files directly
df1 = pd.read_excel('CARD_results.xls')
#remove sheet 2 from Gray et al SI table 2
df2 = pd.read_excel('es0c03803_si_002.xls')

df3 = pd.read_excel('CheckM Results.xlsx')
checkM_dict = mkCheckM(df3)
checkM_block_list = filterCheckM(checkM_dict)

# print(checkM_block_list)




#old phenotype data by class
#df2 = pd.read_excel('phenotype data.xls')

phenotype_dict = {}
rows = []
rows=mkrows(df2)
random.shuffle(rows)

phenotype_dict=mkphenofromrow(rows)
print("before remove phenotype_dict have {} keys".format(len(phenotype_dict)))
for key in checkM_block_list:
    if key in phenotype_dict:
        phenotype_dict.pop(key)
print("after remove phenotype_dict have {} keys".format(len(phenotype_dict)))


#Is this used?
#resistance_dict = {}

CARD_dict = {}
rows = []
CARD_dict=mkCARDdict(df1)
print("before remove CARD_dict have {} keys".format(len(CARD_dict)))
for key in checkM_block_list:
    if key in CARD_dict:
        CARD_dict.pop(key)
print("after remove CARD_dict have {} keys".format(len(CARD_dict)))


#Is this used? If so, where?
import seaborn as sn

df=mkconfusion(CARD_dict,phenotype_dict,pd)
