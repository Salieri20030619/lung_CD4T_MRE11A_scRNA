from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/four_axis_screen_20261017_v15';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3'];rd=lambda f:pd.read_csv(f,sep='\t');a=rd(O/'04_macrophage/consolidated_all_links.tsv.gz');a=a[(a.gene=='DUSP5')&(a.scope=='Treg')&a.primary_pass];a.to_csv(O/'05_intersection/DUSP5_Treg_primary_links.tsv',sep='\t',index=False)
rows=[];complexrows=[]
for lin in ['Mac','CD4']:
 md=rd(M/f'comm_input/{lin}_cells.tsv');gn=rd(M/f'comm_input/{lin}_genes.tsv');g=pd.Index(gn.canonical);x=mmread(M/f'comm_input/{lin}_counts.mtx.gz').T.tocsr()
 def val(gg):
  ix=g.get_indexer(gg);ix=ix[ix>=0];return x[:,ix].toarray() if len(ix) else np.zeros((len(md),0))
 mye=(val(['CSF1R','LYZ2','TYROBP','FCER1G','LST1'])>0).sum(1);epi=(val(['SFTPC','SFTPB','SCGB1A1','EPCAM','KRT8','KRT18'])>0).sum(1);tt=(val(['CD3D','CD3E','CD3G','TRAC'])>0).sum(1)
 query=list(a.ligand.unique()) if lin=='Mac' else ['DUSP5','BMPR1A','BMPR2','TGFBR1','TGFBR2','CDH1','IGF1R','ITGAE','ITGB7','KLRG1','PTPRF','APLP2','CD8A','CD8B1','KLRD1','B2M']
 for sc in (['AMlike','AllMac'] if lin=='Mac' else ['Treg']):
  for s in S:
   mask=md['sample'].eq(s).to_numpy()&(np.ones(len(md),bool) if sc=='AllMac' else md.subtype.eq(sc).to_numpy());n=int(mask.sum())
   for gene in query:
    if gene not in g:continue
    v=val([gene]).ravel();pos=mask&(v>0);rows.append(dict(lineage=lin,scope=sc,sample=s,gene=gene,cells=n,positive=int(pos.sum()),ge2=int((mask&(v>=2)).sum()),UMI=int(v[mask].sum()),CP10k=float(np.mean(v[mask]*1e4/md.loc[mask,'library'])),positive_mye2=int((pos&(mye>=2)).sum()),positive_epi2=int((pos&(epi>=2)).sum()),positive_T2=int((pos&(tt>=2)).sum())))
   if lin=='CD4':
    for name,gg in {'BMP_receptors':['BMPR1A','BMPR2'],'TGF_receptors':['TGFBR1','TGFBR2'],'CD103_integrin':['ITGAE','ITGB7'],'CD8_ab':['CD8A','CD8B1'],'CDH1_KLRG1_receiver':['KLRG1']}.items():
     v=val(gg);ok=v.shape[1]==len(gg);joint=(v>0).all(1) if ok else np.zeros(len(md),bool);complexrows.append(dict(scope=sc,sample=s,complex=name,genes=';'.join(gg),cells=n,measured=ok,joint_positive=int((mask&joint).sum()),pct=100*(mask&joint).sum()/n))
r=pd.DataFrame(rows);r.to_csv(O/'06_novelty/DUSP5_link_identity_counts.tsv',sep='\t',index=False);pd.DataFrame(complexrows).to_csv(O/'06_novelty/DUSP5_receptor_complexes.tsv',sep='\t',index=False)
print('SOURCE_IDENTITY');print(r[r.lineage=='Mac'].to_string(index=False));print('RECEPTORS');print(pd.DataFrame(complexrows).to_string(index=False));print('DUSP5');print(r[r.gene=='DUSP5'].to_string(index=False))