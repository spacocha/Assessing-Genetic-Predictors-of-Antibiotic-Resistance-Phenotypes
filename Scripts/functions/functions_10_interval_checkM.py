#define functions
#ideally, just import these
#however, for now just copy/paste

def mkrows(dfdef):
 rows = []
 for index,row in dfdef.iterrows():
  if index <=2:
    continue
  rows.append((index,row))
 return(rows)

def mkCheckM(df):
  checkM_dict = {}
  for index,row in df.iterrows():
    if index <=1:
      continue
    binId = row["bin id"]
    completeness = row["Completeness"]
    contamination = row["Contamination"]
    checkM_dict[binId] = [completeness,contamination]
  return checkM_dict

def filterCheckM(checkM_dict):
  completenessTH = 90
  contaminationTH = 10
  block_list = []
  for key in checkM_dict:
    values = checkM_dict[key]
    if values[0] <= completenessTH or values[1] >= contaminationTH:
      block_list.append(key.rsplit('_', 1)[0])
  return block_list

def mkphenofromrow(rows):
 phenotype_dict = {}
 for index, row in rows:
  isolateId = row['isolateID']
  phenotype_dict[isolateId] = []
  phenotype_dict[isolateId].append(row['amp_res'])
  phenotype_dict[isolateId].append(row['cip_res'])
  phenotype_dict[isolateId].append(row['tet_res'])
  phenotype_dict[isolateId].append(row['c_res'])
  phenotype_dict[isolateId].append(row['gm_res'])
  phenotype_dict[isolateId].append(row['azm_res'])
  phenotype_dict[isolateId].append(row['cl_res'])
  phenotype_dict[isolateId].append(row['er_res'])
  phenotype_dict[isolateId].append(row['tmp_res'])
  
 return(phenotype_dict)

#This doesn't always work, not sure why
#Use mkrows and mkphenofromrow instead
def mkphenodict(dfdef):
 for index,row in dfdef.iterrows():
  if index <=2:
    continue
  rows.append((index,row))
 for index, row in rows:
  isolateId = row['isolateID']
  phenotype_dict[isolateId] = []
  phenotype_dict[isolateId].append(row['amp_res'])
  phenotype_dict[isolateId].append(row['cip_res'])
  phenotype_dict[isolateId].append(row['tet_res'])
  phenotype_dict[isolateId].append(row['c_res'])
  phenotype_dict[isolateId].append(row['gm_res'])
  phenotype_dict[isolateId].append(row['azm_res'])
  phenotype_dict[isolateId].append(row['cl_res'])
  phenotype_dict[isolateId].append(row['er_res'])
  phenotype_dict[isolateId].append(row['tmp_res'])
 return(phenotype_dict)

def mkCARDdict(dfdef):
 CARD_dict = {}
 rows = []
 for index,row in dfdef.iterrows():
  rows.append((index,row))
 for index, row in rows:
  isolateId = row['isolateID'].rsplit('_', 1)[0]
  if isolateId not in CARD_dict:
    CARD_dict[isolateId] = []
  data = []
  data.append(row['Percentage Length of Reference Sequence'])
  data.append(row['Best_Identities'])
  data.append(row['Antibiotic'])
  CARD_dict[isolateId].append(data)
 return(CARD_dict)

def mkconfusion(CARD_dict,phenotype_dict,pd):
 targets = ["ampicillin","ciprofloxacin","tetracycline","chloramphenicol","gentamicin","azithromycin","colistin","erythromycin","trimethoprim"]
 df=[[-1 for i in range(11)] for j in range(11)]
 for length in range(0,101,10):
  for identity in range(0,101,10):
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
        if drugs=="":
          continue
          # TP
	  # I this this needs to be done at the isolate level, not the gene level
	  #also, I think this should be -1 (resistant), 0 (neither), and 1 (sensitive).
	  #Changed from the original script to TP is has gene and is not sensitive
        if target.lower() in drugs.lower() and (phenotype_dict[isolateId][start_index]==-1 or phenotype_dict[isolateId][start_index]==0):
          m1[0]+=1
          # TN
        elif not target.lower() in drugs.lower() and (phenotype_dict[isolateId][start_index]==1):
          m1[1]+=1
          # FP
        elif target.lower() in drugs.lower() and phenotype_dict[isolateId][start_index]==1:
          m1[2]+=1
          # FN
        elif not target.lower() in drugs.lower() and (phenotype_dict[isolateId][start_index]==-1 or phenotype_dict[isolateId][start_index]==0):
          m1[3]+=1
        start_index+=1
    #print(length,identity,m1[0],m1[1],m1[2],m1[3])
     # TP,TN,FP,FN
    df[length//10][identity//10] = m1
 return(df)

def calcmcc(df,math):
 matrix = [[0 for i in range(11)] for j in range(11)]
 for length in range(0,101,10):
  for identity in range(0,101,10):
    m1 = df[length//10][identity//10]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))==0:
      matrix[length//10][identity//10] = 0
    else:
      matrix[length//10][identity//10] = ((tp*tn)-(fp*fn))/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
 return(matrix)

def calcacc(df):
 matrix = [[0 for i in range(11)] for j in range(11)]
 for length in range(0,101,10):
  for identity in range(0,101,10):
    m1 = df[length//10][identity//10]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tn+fp+fn+tp)==0:
      matrix[length//10][identity//10] = 0
    else:
      matrix[length//10][identity//10] = (tn+tp)/(tn+fp+fn+tp)
 return(matrix)

#Fix the functions so it's one place
#where everything is calculated
#but return just one metric
def calcpre(df):
 matrix = [[0 for i in range(11)] for j in range(11)]
 for length in range(0,101,10):
  for identity in range(0,101,10):
    m1 = df[length//10][identity//10]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tp+fp)==0:
      matrix[length//10][identity//10] = 0
    else:
      matrix[length//10][identity//10] = tp/(tp+fp)
 return(matrix)

def calcrecall(df):
 matrix = [[0 for i in range(11)] for j in range(11)]
 for length in range(0,101,10):
  for identity in range(0,101,10):
    m1 = df[length//10][identity//10]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tp+fn)==0:
      matrix[length//10][identity//10] = 0
    else:
      matrix[length//10][identity//10] = tp/(tp+fn)
 return(matrix)

def calcf1(df):
 matrix = [[0 for i in range(11)] for j in range(11)]
 for length in range(0,101,10):
  for identity in range(0,101,10):
    m1 = df[length//10][identity//10]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tp+0.5*(fp+fn))==0:
      matrix[length//10][identity//10] = 0
    else:
      matrix[length//10][identity//10] = tp/(tp+0.5*(fp+fn))
 return(matrix)

def calcspe(df):
 matrix = [[0 for i in range(11)] for j in range(11)]
 for length in range(0,101,10):
  for identity in range(0,101,10):
    m1 = df[length//10][identity//10]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tn+fp)==0:
      matrix[length//10][identity//10] = 0
    else:
      matrix[length//10][identity//10] = tn/(tn+fp)
 return(matrix)


def mkplot(matrix, metric, pd, plt, sn):
 matrix = pd.DataFrame.from_dict(matrix)
 x_label = []
 y_label = []
 for x in range(0,101,10):
  x_label.append(x)
 for y in range(0,101,10):
  y_label.append(y)
 matrix = pd.DataFrame.from_dict(matrix)
 ax = sn.heatmap(matrix, cmap="YlGnBu")
 ax.invert_yaxis()
 ax.set_xticklabels(x_label)
 ax.set_yticklabels(y_label)
 plt.title(metric)
 plt.xlabel('Identity')
 plt.ylabel('Length')
 #plt.show()
 filename="%s.png" % metric
 plt.savefig(filename,format='png')

def mkplotRedGreen(matrix, metric, pd, plt, sn,current_time):
    # Convert the dictionary to a DataFrame
    matrix = pd.DataFrame.from_dict(matrix)
    x_label = []
    y_label = []
    for x in range(0,101,10):
      x_label.append(x)
    for y in range(0,101,10):
      y_label.append(y)
    # Create binary colormap: 0 → green, 1 → red
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["blue", "red"])

    # Create heatmap
    ax = sn.heatmap(matrix, cmap=cmap, cbar=False, linewidths=0.5, linecolor='gray')

    # Invert y-axis for consistency
    ax.invert_yaxis()
    ax.set_xticklabels(x_label)
    ax.set_yticklabels(y_label)
    # Labeling
    plt.title(metric)
    plt.xlabel('Identity')
    plt.ylabel('Length')

    # Save the figure
    import os
    if not os.path.exists(current_time):
      os.makedirs(current_time)
    filename = f"{current_time}/{metric}.png"
    plt.savefig(filename, format='png')
    plt.close()

def mkbootstrap(obsmat, reps, metric, pd, random, math):
 mcc_maxrmat= [[0 for i in range(11)] for j in range(11)]
 acc_maxrmat= [[0 for i in range(11)] for j in range(11)]
 f1_maxrmat= [[0 for i in range(11)] for j in range(11)]
 pre_maxrmat= [[0 for i in range(11)] for j in range(11)]
 recall_maxrmat= [[0 for i in range(11)] for j in range(11)]
 spe_maxrmat= [[0 for i in range(11)] for j in range(11)]

 
 from datetime import datetime
 for i in range(0,reps):
  print(f"Running {i+1} iteration, {datetime.now()}")
  #Create alpha random matrices in the program
  #Read in real phenotypes
  #remove sheet 2 from Gray et al SI table 2
  df2 = pd.read_excel('es0c03803_si_002.xls')
  #shuffle each of the important antibiotics
  # # use df sample to avoid SettingWithCopyWarning
  # df2['amp_res'] = df2['amp_res'].sample(frac=1).reset_index(drop=True)
  # df2['cip_res'] = df2['cip_res'].sample(frac=1).reset_index(drop=True)
  # df2['tet_res'] = df2['tet_res'].sample(frac=1).reset_index(drop=True)
  # df2['c_res'] = df2['c_res'].sample(frac=1).reset_index(drop=True)
  # df2['gm_res'] = df2['gm_res'].sample(frac=1).reset_index(drop=True)
  # df2['azm_res'] = df2['azm_res'].sample(frac=1).reset_index(drop=True)
  # df2['cl_res'] = df2['cl_res'].sample(frac=1).reset_index(drop=True)

  df1 = pd.read_excel('CARD_results.xls')
  CARD_dict = {}
  rows = []
  df1['Antibiotic'] = df1['Antibiotic'].sample(frac=1).reset_index(drop=True)
  # df1['Percentage Length of Reference Sequence'] = df1['Percentage Length of Reference Sequence'].sample(frac=1).reset_index(drop=True)
  # df1['Best_Identities'] = df1['Best_Identities'].sample(frac=1).reset_index(drop=True)
  CARD_dict=mkCARDdict(df1)

  # random.shuffle(df2.amp_res)
  # random.shuffle(df2.cip_res)
  # random.shuffle(df2.tet_res)
  # random.shuffle(df2.c_res)
  # random.shuffle(df2.gm_res)
  # random.shuffle(df2.azm_res)
  # random.shuffle(df2.cl_res)

  #clear variables
  phenotype_dict = {}
  rows = []
  #remake the phenotype_dict
  rows=mkrows(df2)
  random.shuffle(rows)
  phenotype_dict=mkphenofromrow(rows)
  #make confusion matrix
  df=mkconfusion(CARD_dict,phenotype_dict,pd)
  #calculate one of the metrics from type
  mcc_matrix=calcmcc(df, math)
  acc_matrix=calcacc(df)
  pre_matrix=calcpre(df)
  spe_matrix=calcspe(df)
  f1_matrix=calcf1(df)
  recall_matrix=calcrecall(df)
  
  #max value of random matrices compared to real values
  compareMatrix(mcc_matrix,mcc_maxrmat,obsmat["mcc"])
  compareMatrix(f1_matrix,f1_maxrmat,obsmat["f1"])
  compareMatrix(acc_matrix,acc_maxrmat,obsmat["acc"])
  compareMatrix(pre_matrix,pre_maxrmat,obsmat["pre"])
  compareMatrix(spe_matrix,spe_maxrmat,obsmat["spe"])
  compareMatrix(recall_matrix,recall_maxrmat,obsmat["recall"])
 return {
  "mcc": mcc_maxrmat,
  "f1": f1_maxrmat,
  "acc": acc_maxrmat,
  "pre": pre_maxrmat,
  "spe": spe_maxrmat,
  "recall": recall_maxrmat
}

def compareMatrix(matrix,maxrmat,obsmat):
  for length in range(0,101,10):
    for identity in range(0,101,10):
      if matrix[length//10][identity//10]>=obsmat[length//10][identity//10]:
        maxrmat[length//10][identity//10]+=1
  return(maxrmat)