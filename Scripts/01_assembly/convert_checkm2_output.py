from pathlib import Path
import re
data = []

for path in Path("checkm2").glob('*_*/quality_report.tsv'):

    path_in_str = str(path)
    f = open(path_in_str,"r").readlines()[1:]
    for i in range(len(f)):
      m=re.search('checkm2/(.*)/quality_report.tsv',path_in_str)
      isolateId=m.group(1)
      data.append([isolateId]+f[i].strip().split("\t"))


    

import pandas as pd
df = pd.DataFrame(data, columns=['IsolateId','Name','Completeness','Contamination','Completeness_Model_Used','Translation_Table_Used','Coding_Density','Contig_N50','Average_Gene_Length','Genome_Size','GC_Content','Total_Coding_Sequences','Total_Contigs','Max_Contig_Length','Additional_Notes'])

df.to_csv("checkm2_results.csv")
