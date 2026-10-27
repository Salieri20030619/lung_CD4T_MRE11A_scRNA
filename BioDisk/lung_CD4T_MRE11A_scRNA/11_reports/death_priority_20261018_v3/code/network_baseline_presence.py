from pathlib import Path
import pandas as pd
from pyscenic.cli.utils import load_signatures
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/death_priority_20261018_v2';O=P/'11_reports/death_priority_20261018_v3'
aliases={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(g):return aliases.get(g.upper(),g.upper())
rows=[]
for arm in ['All','BG1']:
 for r in load_signatures(str(V/'joint_network/Mouse'/arm/'seed17/regulons.csv')):
  for g in r.gene2weight:rows.append(dict(arm=arm,tf=canon(r.transcription_factor),target=canon(g),baseline_seed17_present=True))
pd.DataFrame(rows).drop_duplicates().to_csv(O/'tables/network_baseline_seed17_edges.tsv.gz',sep='\t',index=False,encoding='utf-8')
print('BASELINE_SEED17_EDGE_IDENTITY_EXPORTED',len(rows))
