import pandas as pd
data = pd.read_excel('CARD_results.xlsx')
col = ["isolateID", "Best_Hit_ARO", "Best_Identities", "Percentage Length of Reference Sequence", "Drug Class", "Resistance Mechanism", "AMR Gene Family", "Antibiotic"]
value_counts_dict = {}
for i in col:
   value_counts_dict[f"{i}_counts"] = data[i].value_counts().reset_index()

for name, df in value_counts_dict.items():
    df["percentage"] = round(df["count"] / df["count"].sum()*100,2)
    df.to_csv(f"{name}.csv", index=False)