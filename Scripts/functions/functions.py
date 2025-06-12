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
 targets = ["ampicillin","ciprofloxacin","tetracycline","chloramphenicol","gentamicin","azithromycin","colistin"]
 df=[[-1 for i in range(101)] for j in range(121)]
 for length in range(0,121,1):
  for identity in range(0,101,1):
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
	  # I this this needs to be done at the isolate level, not the gene level
	  #also, I think this should be -1 (resistant), 0 (neither), and 1 (sensitive).
	  #Changed from the original script to TP is has gene and is not sensitive
        if target.lower() in drugs.lower() and (phenotype_dict[isolateId][start_index]==-1 or phenotype_dict[isolateId][start_index]==0):
          m1[0]+=1
          # TN
        if not target.lower() in drugs.lower() and (phenotype_dict[isolateId][start_index]==1):
          m1[1]+=1
          # FP
        if target.lower() in drugs.lower() and phenotype_dict[isolateId][start_index]==1:
          m1[2]+=1
          # FN
        if not target.lower() in drugs.lower() and (phenotype_dict[isolateId][start_index]==-1 or phenotype_dict[isolateId][start_index]==0):
          m1[3]+=1
        start_index+=1
    #print(length,identity,m1[0],m1[1],m1[2],m1[3])
     # TP,TN,FP,FN
    df[length][identity] = m1
 return(df)

def calcmcc(df,math):
 matrix = [[0 for i in range(101)] for j in range(121)]
 for length in range(0,121,1):
  for identity in range(0,101,1):
    m1 = df[length][identity]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))==0:
      matrix[length][identity] = 0
    else:
      matrix[length][identity] = ((tp*tn)-(fp*fn))/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
 return(matrix)

def calcacc(df):
 matrix = [[0 for i in range(101)] for j in range(121)]
 for length in range(0,121,1):
  for identity in range(0,101,1):
    m1 = df[length][identity]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tn+fp+fn+tp)==0:
      matrix[length][identity] = 0
    else:
      matrix[length][identity] = (tn+tp)/(tn+fp+fn+tp)
 return(matrix)

#Fix the functions so it's one place
#where everything is calculated
#but return just one metric
def calcpre(df):
 matrix = [[0 for i in range(101)] for j in range(121)]
 for length in range(0,121,1):
  for identity in range(0,101,1):
    m1 = df[length][identity]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tp+fp)==0:
      matrix[length][identity] = 0
    else:
      matrix[length][identity] = tp/(tp+fp)
 return(matrix)

def calcrecall(df):
 matrix = [[0 for i in range(101)] for j in range(121)]
 for length in range(0,121,1):
  for identity in range(0,101,1):
    m1 = df[length][identity]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tp+fn)==0:
      matrix[length][identity] = 0
    else:
      matrix[length][identity] = tp/(tp+fn)
 return(matrix)

def calcf1(df):
 matrix = [[0 for i in range(101)] for j in range(121)]
 for length in range(0,121,1):
  for identity in range(0,101,1):
    m1 = df[length][identity]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tp+0.5*(fp+fn))==0:
      matrix[length][identity] = 0
    else:
      matrix[length][identity] = tp/(tp+0.5*(fp+fn))
 return(matrix)

def calcspe(df):
 matrix = [[0 for i in range(101)] for j in range(121)]
 for length in range(0,121,1):
  for identity in range(0,101,1):
    m1 = df[length][identity]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tn+fp)==0:
      matrix[length][identity] = 0
    else:
      matrix[length][identity] = tn/(tn+fp)
 return(matrix)


def mkplot(matrix, metric, pd, plt, sn):
 matrix = pd.DataFrame.from_dict(matrix)
 x_label = []
 y_label = []
 for x in range(101):
  x_label.append(x)
 for y in range(121):
  y_label.append(y)
 matrix = pd.DataFrame.from_dict(matrix)
 ax = sn.heatmap(matrix, cmap="YlGnBu")
 ax.invert_yaxis()
 plt.title(metric)
 plt.xlabel('Identity')
 plt.ylabel('Length')
 #plt.show()
 filename="%s.png" % metric
 plt.savefig(filename,format='png')

def mkbootstrap(CARD_dict, obsmat, reps, metric, pd, random, math):
 maxrmat= [[0 for i in range(101)] for j in range(121)]
 for i in range(0,reps):
  #Create alpha random matrices in the program
  #Read in real phenotypes
  #remove sheet 2 from Gray et al SI table 2
  df2 = pd.read_excel('es0c03803_si_002.xls')
  #shuffle each of the important antibiotics
  random.shuffle(df2.amp_res)
  random.shuffle(df2.cip_res)
  random.shuffle(df2.tet_res)
  random.shuffle(df2.c_res)
  random.shuffle(df2.gm_res)
  random.shuffle(df2.azm_res)
  random.shuffle(df2.cl_res)
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
  if metric=="MCC":
   matrix=calcmcc(df, math)
  elif metric=="ACC":
   matrix=calcacc(df)
  elif metric=="PRE":
   matrix=calcpre(df)
  elif metric=="SPE":
   matrix=calcspe(df)
  elif metric=="F1":
   matrix=calcf1(df)
  elif metric=="RECALL":
   matrix=calcrecall(df)
  
  #max value of random matrices compared to real values
  for length in range(0,121,1):
   for identity in range(0,101,1):
    if matrix[length][identity]>=obsmat[length][identity]:
     maxrmat[length][identity]=maxrmat[length][identity]+1
  return(maxrmat)


