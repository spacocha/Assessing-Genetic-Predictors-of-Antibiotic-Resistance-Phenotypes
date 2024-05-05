# find gene_identifier to isolate_id
f =  open("all_16s.fasta","r")
pair = {}
for line in f:
 if ">" in line:
     gene_identifier = line.split(" ")[0][1:]
     isolate_ID = line.split(" ")[1]
     pair[gene_identifier] = isolate_ID

# find isolate_id to srr
import pandas as pd
df = pd.read_excel("data.xlsx")
pair_excel = {}
iso_to_srr = {}
for row in df.itertuples():
 isolate_ID = row[1]
 srr = row[2]
 pair_excel[srr] = isolate_ID
 iso_to_srr[isolate_ID] = srr

# find gene_identifier to srr
result_pair = {}
for gene_identifier in pair.keys():
 isolate_ID = pair[gene_identifier]
 isolate_ID_short = "_".join(isolate_ID.split('_')[:-1])
 find_result = pair_excel.get(isolate_ID_short,None)
 if not find_result:
   find_result = pair_excel.get(isolate_ID,None)
 if  find_result:
   result_pair[gene_identifier] = find_result
print(result_pair)
