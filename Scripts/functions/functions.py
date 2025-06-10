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

#Fix the functions so it's one place
#where everything is calculated
#but return just one metric
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
      pre = tp/(tp+fp)
      matrix[length][identity] = pre
 return(matrix)


def mkplot(matrix, type, pd, plt, sn):
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
 plt.title(type)
 plt.xlabel('Identity')
 plt.ylabel('Length')
 #plt.show()
 filename="%s.png" % type
 plt.savefig(filename,format='png')


