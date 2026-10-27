from pathlib import Path
import pandas as pd,h5py,json
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/integrated_macrophage_CD4_20261015_v4'
for f in (P/'02_analysis/scenic2_20261001_092655_ZKNWD6/tables').iterdir():print('SCENIC_TABLE',f.name)
for f in (P/'02_analysis/scenic1_20260929_103513_xRoEnr/Mouse/All').glob('**/*'):
 if f.is_file() and ('regulon' in f.name or 'signature' in f.name):print('REG',f)
a=pd.read_csv(V/'03_receiver/TF_all_regulon_sample_activity.tsv.gz',sep='\t');print('PROX_TFS',sorted(a.loc[a.tf.isin(['STAT5A','STAT5B','NFKB1','NFKB2','RELA','RELB','JUN','JUNB','FOS','FOSB','STAT1','IKZF2']), 'regulon'].unique()))
for rel in ['stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad','stage4_20260905_215942_3XPis4/objects/TNK_harmony.h5ad']:
 with h5py.File(P/'02_analysis/GSE314943'/rel) as h:
  for key in ['Broad','State','H06']:
   if key in h['obs']:
    x=h['obs'][key];print(rel,key,list(x['categories'].asstr()[:]) if isinstance(x,h5py.Group) else list(set(x.asstr()[:])))
for f in [P/'11_reports/methods_complete_20261016_v13/supplement/S04_death_program_comparisons.tsv',P/'11_reports/methods_complete_20261016_v13/supplement/S05_mito_autophagy_comparisons.tsv']:
 a=pd.read_csv(f,sep='\t');print(f.name,list(a),a.head(2).to_dict('records'))