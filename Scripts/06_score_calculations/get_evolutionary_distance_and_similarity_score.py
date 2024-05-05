from ete3 import PhyloTree
import pandas as pd
import math


# read the score for pair of isolates
df = pd.read_excel('data.xlsx',0)
results = {}
results_str = {}
resistances = ["amp_res","amx_res","ctx_res","er_res","azm_res","tet_res","gm_res","cip_res","cc_res","tmp_res","nitro_res","cl_res","c_res"]
for i in range(len(df['SRR number'].tolist())):
 srr = df['SRR number'].tolist()[i]
 if str(srr) == 'nan':
   continue
 results[srr] = []
 results_str[srr]=""
 for v in resistances:
   num = df[v].tolist()[i]
   results[srr].append(num)
   if num >=0:
     results_str[srr]+="1"
   elif num == -1:
     results_str[srr]+="0"
   else:
     results_str[srr]+="n"
print(results_str)


from ete3 import PhyloTree
import pandas as pd
import math

def get_score(n1,n2):
 score = 0
 counted = 0
 for i in range(13):
   s1 = results_str[n1][i]
   s2 = results_str[n2][i]
   if s1==s2 and s1!="n":
     score+=1
   if s1!="n" or s2!="n":
     counted+=1
 return score/counted

tree = PhyloTree("phyml_tree.txt")

done_pair = []
# Get the two species (nodes) from the tree
distance_dict=[]
score_dict=[]
for key1 in result_pair.keys():
 short_key1 = key1[:10]
 key1_name = result_pair[key1]
 gene1 = tree&short_key1
 done = []

 for key2 in result_pair:
   short_key2 = key2[:10]
   if key1+key2 in done_pair or key2+key1 in done_pair:
     continue
   done_pair.append(key1+key2)
   if short_key2 in done:
     continue
   done.append(short_key2)
   if short_key2 == short_key1:
     continue

   gene2 = tree&short_key2
   distance = round(gene1.get_distance(gene2),5)

   # Print the distance between them
   distance_dict.append(str(distance))
   score_dict.append(str(get_score(key1_name,d[key2])))
print(distance_dict)
print(score_dict)
df = pd.DataFrame({
   "distance":distance_dict,
   "count":score_dict
})
display(df)
