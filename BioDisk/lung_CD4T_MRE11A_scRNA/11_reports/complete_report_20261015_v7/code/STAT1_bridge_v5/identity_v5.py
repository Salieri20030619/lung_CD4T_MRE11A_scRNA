from pathlib import Path
import sys,json
import pandas as pd,numpy as np
from scipy.io import mmread
O=Path(sys.argv[1]);P=O.parents[1];V=O.parent/'integrated_macrophage_CD4_20261015_v4';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';ns={};exec((V/'code/mouse_stats_v4.py').read_text().split('phen=rd(')[0],ns);S=ns['SAMPLES'];assoc=ns['assoc'];wr=lambda a,n:a.to_csv(O/'tables'/n,sep='\t',index=False,na_rep='NA');G=json.loads((O/'provenance/analysis_spec.json').read_text())['genes'];rows=[];results=[];source=[]
for lin in ['CD4','Mac']:
 md=pd.read_csv(M/'comm_input'/f'{lin}_cells.tsv',sep='\t');gn=pd.read_csv(M/'comm_input'/f'{lin}_genes.tsv',sep='\t');gi=pd.Index(gn.canonical);x=mmread(M/'comm_input'/f'{lin}_counts.mtx.gz').T.tocsr()
 def detected(gs):
  ix=gi.get_indexer(gs);ix=ix[ix>=0];return np.asarray((x[:,ix]>0).sum(axis=1)).ravel()
 my=detected(['CSF1R','LST1','TYROBP','FCER1G','LYZ2'])>=2;ep=detected(['EPCAM','SFTPA1','SFTPB','SFTPC','SFTPD','SCGB1A1'])>=2;tt=detected(['CD3D','CD3E','TRAC','CD247'])>=2
 if lin=='CD4':
  auc=pd.read_csv(P/'02_analysis/scenic2_20261001_092655_ZKNWD6/objects/Mouse_All_17_DropOverlap_auc.csv.gz',index_col=0).loc[md.cell,'Stat1(+)'].to_numpy()
  for scope in ['Core','Naive','ActCD4','Treg']:
   xx=[];yy=[]
   for sample in S:
    base=md['sample'].eq(sample).to_numpy()&(np.ones(len(md),bool) if scope=='Core' else md.subtype.eq(scope).to_numpy());keep=base&~my&~ep;xx.append(auc[keep].mean() if keep.any() else np.nan);gvals=[]
    for g in G:
     j=gi.get_indexer([g])[0];vals=x[:,j].toarray().ravel() if j>=0 else np.zeros(len(md));positive=base&(vals>0);rows.append(dict(scope=scope,sample=sample,gene=g,cells=int(base.sum()),kept_cells=int(keep.sum()),positive_cells=int(positive.sum()),myeloid_flag_positive=int((positive&my).sum()),epithelial_flag_positive=int((positive&ep).sum()),T_supported_positive=int((positive&tt).sum())));gvals.append(float(np.mean(vals[keep]*1e4/md.loc[keep,'library'])) if j>=0 and keep.any() else np.nan)
    yy.append(gvals)
   X=pd.DataFrame({'Stat1(+)':xx},index=S);Y=pd.DataFrame(yy,index=S,columns=G).dropna(axis=1)
   if X.notna().all().all():results.append(assoc(X,Y,scope,'flag_exclusion',17,'STAT1_gene_QC'))
 else:
  for scope in ['AMlike','AllMac']:
   for sample in S:
    mask=md['sample'].eq(sample).to_numpy()&(np.ones(len(md),bool) if scope=='AllMac' else md.subtype.eq(scope).to_numpy())
    for name,gs in {'IL12_heterodimer':['IL12A','IL12B'],'IL27_heterodimer':['IL27','EBI3'],'IFNG':['IFNG'],'IFNB1':['IFNB1'],'IFNA1':['IFNA1'],'IFNAB':['IFNAB'],'C3':['C3'],'CXCL16':['CXCL16']}.items():
     ii=gi.get_indexer(gs);ok=(ii>=0).all();positive=mask&(detected(gs)==len(gs)) if ok else np.zeros(len(md),bool);source.append(dict(scope=scope,sample=sample,complex=name,genes=';'.join(gs),measured=ok,cells=int(mask.sum()),positive_cells=int(positive.sum()),pct=100*positive.sum()/mask.sum(),T_flag_positive=int((positive&tt).sum()),epithelial_flag_positive=int((positive&ep).sum()),myeloid_supported_positive=int((positive&my).sum())))
wr(pd.DataFrame(rows),'CD4_executor_identity_audit.tsv');res=pd.concat(results) if results else pd.DataFrame(columns=['scope','model','y','rho','P','BH']);wr(res,'flag_exclusion_STAT1_gene.tsv');wr(pd.DataFrame(source),'source_complex_identity_audit.tsv');print('QC_GSDMD');print(res[(res.model=='Pooled')&(res.y=='GSDMD')][['scope','rho','P','BH']].to_string(index=False));print('SOURCE_CASES');s=pd.DataFrame(source);print(s[s['sample'].str.startswith('Silica')].groupby(['scope','complex']).positive_cells.sum().to_string());print('RETAINED');print(pd.DataFrame(rows).drop_duplicates(['scope','sample'])[['scope','sample','cells','kept_cells']].to_string(index=False))