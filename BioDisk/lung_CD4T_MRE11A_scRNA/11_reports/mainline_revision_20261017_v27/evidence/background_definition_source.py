from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json,scipy.sparse as sp
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/three_cohort_four_axis_20261017_v16';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/comm_input'
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(x):
 x=str(x).strip().upper();return alias.get(x,x)
inv=[];mapping=[]
def project(raw,sy):
 sy=np.array([canon(s) for s in sy]);valid=(sy!='')&~np.char.startswith(sy,'ENS');genes=np.array(sorted(set(sy[valid])));pr=sp.csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),pd.Index(genes).get_indexer(sy[valid]))),shape=(raw.shape[1],len(genes)));return (raw@pr).tocsr(),genes,sy

def export(ds,scope,x,g,m):
 assert x.shape==(len(m),len(g));assert (x.data>=0).all() and np.allclose(x.data,np.rint(x.data))
 samples=sorted(m['sample'].unique());counts=[];pct=[];met=[]
 for s in samples:
  take=m['sample'].eq(s).to_numpy();z=x[take];counts.append(np.asarray(z.sum(0)).ravel());pct.append(np.asarray((z>0).sum(0)).ravel()/sum(take));a=m.loc[take].iloc[0];met.append(dict(sample=s,dataset=ds,scope=scope,group=a['group'],day=a['day'],cells=int(sum(take)),library=float(m.loc[take,'library'].sum())))
 stem=ds+'__'+scope
 pd.DataFrame(np.array(counts).T,index=pd.Index(g,name='gene'),columns=samples).astype('int64').to_csv(O/'inputs'/f'{stem}.counts.tsv.gz',sep='\t')
 pd.DataFrame(np.array(pct).T,index=pd.Index(g,name='gene'),columns=samples).to_csv(O/'inputs'/f'{stem}.fraction.tsv.gz',sep='\t')
 pd.DataFrame(met).to_csv(O/'inputs'/f'{stem}.samples.tsv',sep='\t',index=False)
 inv.extend(met);print('EXPORTED',stem,len(m),len(g),len(samples),flush=True)
 # Record cell-level background and selected lineage signals for later sensitivity, not hard exclusion.
 epi={'SFTPA1','SFTPA2','SFTPB','SFTPC','SFTPD','SCGB1A1','SCGB3A1','SCGB3A2','ABCA3','SLC34A2','NAPSA'};mye={'LYZ','LYZ1','LYZ2','S100A8','S100A9','CSF1R','LST1','TYROBP','FCER1G'}
 import re
 bg=np.array([v in epi|mye or bool(re.match(r'^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',v)) for v in g]);m=m.copy();m['background_fraction']=np.asarray(x[:,bg].sum(1)).ravel()/m.library.to_numpy();m.to_csv(O/'inputs'/f'{stem}.cells.tsv.gz',sep='\t',index=False)
 if scope in ['Core','TLike']:
  sp.save_npz(O/'inputs'/f'{stem}.cellcounts.npz',x);pd.Series(g).to_csv(O/'inputs'/f'{stem}.genes.txt',index=False,header=False)
for lin in ['CD4','Mac']:
 m=pd.read_csv(M/f'{lin}_cells.tsv',sep='\t');gs=pd.read_csv(M/f'{lin}_genes.tsv',sep='\t');raw=mmread(M/f'{lin}_counts.mtx.gz').T.tocsr();x,g,sy=project(raw,gs.canonical);m['day']=140
 for orig,ca in zip(gs.canonical,sy):mapping.append(dict(dataset='GSE314943',lineage=lin,source=orig,canonical=ca))
 scopes=['Core','Naive','ActCD4','Treg'] if lin=='CD4' else ['AllMac','AMlike']
 print('SUBTYPES',lin,m.subtype.value_counts().to_dict(),flush=True)
 for scope in scopes:
  take=np.ones(len(m),bool) if scope in ['Core','AllMac'] else m.subtype.eq(scope).to_numpy()
  export('GSE314943',scope,x[take],g,m.loc[take].reset_index(drop=True))
for ds,path in [('GSE183682','02_analysis/screen2_CD4_target_wz1lfx/GSE183682/CD4_core.h5ad'),('GSE274372','02_analysis/GSE274372/stage9_20260908_140721_Yj56oq/objects/CD4_author_reclustered.h5ad'),('GSE274372Mac','02_analysis/GSE274372/stage3_20260907_171359_MgRRUd/objects/whole_annotated.h5ad')]:
 a=ad.read_h5ad(P/path);raw=sp.csr_matrix(a.layers['counts']);sy=a.var.symbol.astype(str).values if 'symbol' in a.var else a.var_names.to_numpy();x,g,sy=project(raw,sy)
 for orig,ca in zip(a.var_names,sy):mapping.append(dict(dataset=ds,lineage='Mac' if ds.endswith('Mac') else 'CD4',source=orig,canonical=ca))
 m=pd.DataFrame({'cell':a.obs_names,'sample':a.obs['sample'].astype(str).to_numpy(),'library':np.asarray(raw.sum(1)).ravel()})
 if ds=='GSE183682':
  mp={'GSM5567501':('Saline',7),'GSM5567502':('Silica',7),'GSM5567503':('Saline',56),'GSM5567504':('Silica',56)};m['group']=m['sample'].map({k:v[0] for k,v in mp.items()});m['day']=m['sample'].map({k:v[1] for k,v in mp.items()});m['subtype']=a.obs.core_label.astype(str).replace({'NaiveCD4':'Naive'}).to_numpy()
  for scope in ['Core','Naive','ActCD4','Treg']:
   take=np.ones(len(m),bool) if scope=='Core' else m.subtype.eq(scope).to_numpy();export(ds,scope,x[take],g,m.loc[take].reset_index(drop=True))
 else:
  m['day']=pd.to_numeric(a.obs.day).to_numpy();m['group']='D'+m.day.astype(int).astype(str)
  if ds.endswith('Mac'):take=a.obs.Broad.astype(str).eq('Mac').to_numpy();scope='AllMac'
  else:take=a.obs.R08.astype(str).isin(['0','1']).to_numpy();scope='TLike';assert take.sum()==494
  m['subtype']=scope;export('GSE274372',scope,x[take],g,m.loc[take].reset_index(drop=True))
pd.DataFrame(inv).to_csv(O/'tables/sample_inventory.tsv',sep='\t',index=False);pd.DataFrame(mapping).to_csv(O/'tables/feature_mapping.tsv.gz',sep='\t',index=False)
(O/'logs/export_complete.json').write_text(json.dumps({'success':True,'sample_scope_rows':len(inv)}),encoding='utf-8')