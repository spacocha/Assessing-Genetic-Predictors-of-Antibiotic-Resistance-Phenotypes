import pandas as pd

df1 = pd.read_excel('CARD_results.xls')
df2=pd.read_excel('es0c03803_si_002.xls')
tcs = []
bc = []
for line in df2.itertuples():
  if line.media=='tcs':
    tcs.append(line.isolateID)
  elif line.media == 'bc':
    bc.append(line.isolateID)
total_bc = 0
total_tcs = 0
have_bc = set()
have_tcs = set()

for line in df1.itertuples():
  isolateId = line.isolateID.rsplit('_', 1)[0]
  if isolateId in bc:
    if 'benzalkonium chloride' in line.Antibiotic.lower():
      have_bc.add(isolateId)
  if isolateId in tcs:
    if 'triclosan' in line.Antibiotic.lower():
        have_tcs.add(isolateId)
    
print(f"total bc: {len(bc)}")
print(f"have bc: {len(have_bc)}")
print(f"total tcs: {len(tcs)}")
print(f"have tcs: {len(have_tcs)}")