# continue from 08_confusion_matrix

import pandas as pd
import math
import random
for i in range(0,100):
  xls = pd.ExcelFile("random/random"+str(i)+".xlsx")
  df1 = pd.read_excel(xls, 'CARD Full results')
  df2 = pd.read_excel(xls, 'phenotype data')
  phenotype_dict = {}
  rows = []
  for index,row in df2.iterrows():
    rows.append((index,row))
  random.shuffle(rows)
  for index, row in rows:
    isolateId = row['isolateID']
    phenotype_dict[isolateId] = []
    phenotype_dict[isolateId].append(row['penam'])
    if (row['macrolide1']+row['macrolide2']) == -2:

      phenotype_dict[isolateId].append(-1)
      phenotype_dict[isolateId].append(-1)
    else:
      phenotype_dict[isolateId].append(1)
      phenotype_dict[isolateId].append(1)
    phenotype_dict[isolateId].append(row['tetracycline'])
    phenotype_dict[isolateId].append(row['aminoglycoside'])
    phenotype_dict[isolateId].append(row['fluoroquinolone'])
    phenotype_dict[isolateId].append(row['sulfonamide'])
    phenotype_dict[isolateId].append(row['peptide'])
    phenotype_dict[isolateId].append(row['phenicol'])
  resistance_dict = {}

  result = []
  groups = {}
  result2  = []
  CARD_dict = {}
  rows = []
  for index,row in df1.iterrows():
    rows.append((index,row))
  random.shuffle(rows)
  for index, row in rows:
    isolateId = row['isolateID'].rsplit('_', 1)[0]
    if isolateId not in CARD_dict:
      CARD_dict[isolateId] = []
    data = []
    data.append(row['Percentage Length of Reference Sequence'])
    data.append(row['Best_Identities'])
    data.append(row['Drug Class'])
    CARD_dict[isolateId].append(data)
  import seaborn as sn

  targets = ["penam","macrolide","macrolide","tetracycline","aminoglycoside","fluoroquinolone","sulfonamide","peptide","phenicol"]
  df=[[-1 for i in range(101)] for j in range(101)]
  for length in range(0,101,1):
    for identity in range(60,101,1):
      m1 = [0,0,0,0]
      for isolateId in CARD_dict:
        drugs = ""
        if isolateId not in phenotype_dict:
          continue
        for data in CARD_dict[isolateId]:
          if data[0] >= length  and data[1] >= identity:
            drugs = drugs + data[2]
        start_index = 0
        for target in targets:
          if pd.isna(drugs):
            continue
            # TP
          if target.lower() in drugs.lower() and (phenotype_dict[isolateId][start_index]==1 or phenotype_dict[isolateId][start_index]==0):
            m1[0]+=1
            # TN
          if not target.lower() in drugs.lower() and (phenotype_dict[isolateId][start_index]==-1):
            m1[1]+=1
            # FP
          if target.lower() in drugs.lower() and phenotype_dict[isolateId][start_index]==-1:
            m1[2]+=1
            # FN
          if not target.lower() in drugs.lower() and (phenotype_dict[isolateId][start_index]==1 or phenotype_dict[isolateId][start_index]==0):
            m1[3]+=1
          start_index+=1
      print(length,identity,m1[0],m1[1],m1[2],m1[3])
      # TP,TN,FP,FN
      df[length][identity] = m1
  with open("random_result/random"+str(i)+".xlsx","w") as f:
    f.write(str(df))



import math
target_matrix = eval(open("target.xlsx","r").read())

import seaborn as sn
import math
import pandas as pd
target_matrix = eval(open("target.xlsx","r").read())
pd.DataFrame(target_matrix)
result_matrix = [[0 for i in range(101)] for j in range(101)]
f1_matrix = [[0 for i in range(101)] for j in range(101)]
for i in range(100):
  with open("random_result/random"+str(i)+".xlsx","r") as f:
    df = eval(f.read())
  for length in range(0,101,1):
    for identity in range(60,101,1):

      m1 = df[length][identity]
      tp = m1[0]
      tn = m1[1]
      fp = m1[2]
      fn = m1[3]
      m2 = target_matrix[length][identity]
      tp2 = m2[0]
      tn2 = m2[1]
      fp2 = m2[2]
      fn2 = m2[3]

      if (tp+fp)==0:
        f1_matrix[length][identity] = 0
      else:
        acc = (tn+tp)/(tn+fp+fn+tp)
        sen = tp/(tp+fn)
        spe = tn/(tn+fp)
        recall = tp/(tp+fn)
        pre = tp/(tp+fp)
        mcc = ((tp*tn)-(fp*fn))/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        f1 = tp/(tp+0.5*(fp+fn))
        acc_o = (tn2+tp2)/(tn2+fp2+fn2+tp2)
        sen_o = tp2/(tp2+fn2)
        spe_o = tn2/(tn2+fp2)
        recall_o = tp2/(tp2+fn2)
        pre_o = tp2/(tp2+fp2)
        mcc_o = ((tp2*tn2)-(fp2*fn2))/math.sqrt((tp2+fp2)*(tp2+fn2)*(tn2+fp2)*(tn2+fn2))
        f1_o = tp2/(tp2+0.5*(fp2+fn2))
        # print(mcc,mcc_o,m1,m2)
        if mcc > mcc_o:
          result_matrix[length][identity] +=1
for length in range(0,101,1):
  for identity in range(60,101,1):
    if result_matrix[length][identity] >=1:
      result_matrix[length][identity] =0
    else:
      result_matrix[length][identity] = 1

import matplotlib.pyplot as plt # data visualization
result_matrix = pd.DataFrame.from_dict(result_matrix)
x_label = []
y_label = []
for x in range(60,101):
  x_label.append(x)
for y in range(101):
  y_label.append(y)
print(result_matrix)
result_matrix = pd.DataFrame(result_matrix,columns=x_label)
ax = sn.heatmap(result_matrix, cmap="YlGn",vmin=0,vmax=1)
ax.invert_yaxis()
plt.title("100 Random Sample 1% P-value bootstrap MCC")
plt.xlabel('Identity')
plt.ylabel('Length')

plt.show()
