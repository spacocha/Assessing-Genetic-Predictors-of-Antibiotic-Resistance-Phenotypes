from pathlib import Path
data = []

for path in Path("rgi_output").glob('*_*/*.txt'):

    path_in_str = str(path)
    f = open(path_in_str,"r").readlines()[1:]
    for i in range(len(f)):
      isolateId = path_in_str.split("/")[-1].split(".")[0].replace("_rgi_output","")
      data.append([isolateId]+f[i].strip().split("\t"))


    

import pandas as pd
df = pd.DataFrame(data, columns=['isolateID','ORF_ID','Contig','Start','Stop','Orientation','Cut_Off','Pass_Bitscore','Best_Hit_Bitscore','Best_Hit_ARO','Best_Identities','ARO','Model_type','SNPs_in_Best_Hit_ARO','Other_SNPs','Drug Class','Resistance Mechanism','AMR Gene Family','Predicted_DNA','Predicted_Protein','CARD_Protein_Sequence','Percentage Length of Reference Sequence','ID','Model_ID','Nudged','Note','Hit_Start','Hit_End','Antibiotic'])

df.to_csv("CARD_results.csv")
