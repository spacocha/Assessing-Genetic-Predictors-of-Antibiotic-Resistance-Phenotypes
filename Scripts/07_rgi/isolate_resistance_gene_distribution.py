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

#If I change the functions confusion_matrix_10_interval_checkM_def
#change this too
checkM_dict = mkCheckM(df3)
checkM_block_list = filterCheckM(checkM_dict)
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

CARD_dict = {}
rows = []
CARD_dict=mkCARDdict(df1)

print("before remove CARD_dict have {} keys".format(len(CARD_dict)))
for key in checkM_block_list:
    if key in CARD_dict:
        CARD_dict.pop(key)
print("after remove CARD_dict have {} keys".format(len(CARD_dict)))
gene_res_dict={}
gene_sens_dict={}
targets = ["trimethoprim","ampicillin","ciprofloxacin","tetracycline","chloramphenicol","gentamicin","azithromycin","colistin","erythromycin"]
for target in targets:
    gene_res_dict[target]=0
    gene_sens_dict[target]=0
for isolateId in phenotype_dict:
    drugs = ""
    if isolateId in CARD_dict:
     for data in CARD_dict[isolateId]:
        drugs = drugs + data[2]
    start_index = 0
    for target in targets:
        if target.lower() in drugs.lower():
            #It has a resistance gene
            gene_res_dict[target]+=1
        else:
            gene_sens_dict[target]+=1

with open('gene_res_dict.txt','w') as file:
 file.write(str(gene_res_dict))

with open('gene_sens_dict.txt','w') as file:
 file.write(str(gene_sens_dict))

