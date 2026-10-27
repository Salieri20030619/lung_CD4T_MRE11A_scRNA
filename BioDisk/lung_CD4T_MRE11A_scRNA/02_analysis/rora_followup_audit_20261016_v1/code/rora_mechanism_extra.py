from pathlib import Path
import sys,hashlib,json
import pandas as pd,numpy as np
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=Path(sys.argv[1]);M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/comm_input';I=P/'02_analysis/GSE274372/stage10_20260908_141459_a8IPaq/tables'
g=pd.read_csv(M/'CD4_genes.tsv',sep='\t');m=pd.read_csv(M/'CD4_cells.tsv',sep='\t');x=mmread(M/'CD4_counts.mtx.gz').T.tocsr();rows=[]
for target,aliases in {'TP53':['TP53','TRP53'],'PRKN':['PRKN','PARK2'],'DRAM1':['DRAM1'],'CAV3':['CAV3']}.items():
 idx=np.flatnonzero(g.canonical.isin(aliases));assert len(idx)<=1
 for scope in ['Core','Naive','ActCD4','Treg']:
  for s in m['sample'].unique():
   take=(m['sample']==s)&((m.subtype==scope) if scope!='Core' else True);ids=np.flatnonzero(take);a=x[ids][:,idx].toarray().ravel() if len(idx) else None
   rows.append(dict(target=target,matched_canonical=';'.join(g.iloc[idx].canonical),matched_mouse_symbol=';'.join(g.iloc[idx].mouse_symbol),scope=scope,sample=s,measured=bool(len(idx)),cells=len(ids),positive_cells=int((a>0).sum()) if a is not None else None,fraction_positive=float((a>0).mean()) if a is not None else None))
pd.DataFrame(rows).to_csv(O/'mechanism_alias_checked.tsv',sep='\t',index=False,encoding='utf-8')
f=pd.read_csv(I/'features.tsv.gz',sep='\t',index_col=0);ct=pd.read_csv(I/'TLike_counts.tsv.gz',sep='\t',index_col=0);meta=pd.read_csv(I/'sample_summary.tsv',sep='\t');meta=meta[meta.scope.eq('TLike')].copy();ix=f.index[f.symbol.str.upper().eq('RORA')];assert len(ix)==1;meta['RORA_UMI']=ct.loc[ix[0],meta['sample']].to_numpy();meta['RORA_raw_CPM']=meta.RORA_UMI/meta.library_UMI*1e6;meta.to_csv(O/'time_RORA_raw_counts.tsv',sep='\t',index=False,encoding='utf-8')
manifest=[]
for n in ['features.tsv.gz','TLike_counts.tsv.gz','sample_summary.tsv']:
 p=I/n;manifest.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
pd.DataFrame(manifest).to_csv(O/'time_input_manifest.tsv',sep='\t',index=False,encoding='utf-8');print(meta[['sample','day','cells','RORA_UMI','RORA_raw_CPM']].to_string(index=False));print('EXTRA_COMPLETE')
