import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json
from scipy import sparse,io
import statslib as S
O=S.O;P=O.parent.parent;B=P/'11_reports/three_cohort_four_axis_20261017_v16';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF/comm_input'
def write(x,n):x.to_csv(O/'tables'/n,sep='\t',index=False)
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
def canon(v):return alias.get(str(v).upper(),str(v).upper())
lr=pd.read_csv(B/'resources/CellChat_protein_LR_expanded.tsv',sep='\t').drop_duplicates('interaction_name').fillna('')
cand=pd.read_csv(O/'tables/frozen_candidate_pairs.tsv',sep='\t');needed=sorted(set(';'.join(lr.ligand_parts).split(';'))|set(';'.join(lr.receptor_parts).split(';'))|set(cand.gene));needed=[g for g in needed if g]
inputs={}
m=pd.read_csv(M/'Mac_cells.tsv',sep='\t');g=pd.read_csv(M/'Mac_genes.tsv',sep='\t').canonical.map(canon).to_numpy();x=io.mmread(M/'Mac_counts.mtx.gz').T.tocsr();m['day']=140;m['dataset']='GSE314943'
inputs['GSE314943']=(x,m,g)
a=ad.read_h5ad(O/'objects/GSE183682_whole_clustered.h5ad');ann=json.loads((O/'inputs/GSE183682_annotation.json').read_text());a.obs['Broad']=a.obs.cluster.astype(str).map(ann);take=a.obs.Broad.eq('Mac').to_numpy();m=a.obs.loc[take].copy().reset_index(names='cell');m['dataset']='GSE183682';inputs['GSE183682']=(a.layers['counts'][take].tocsr(),m,np.array([canon(g) for g in a.var_names]))
rows=[];coverage=[];allcounts=[];allmet=[];strictgenes=['C1QA','C1QB','C1QC','ADGRE1','MERTK','MSR1','MARCO']
for ds,(x,m,g) in inputs.items():
 for variant in ['AllMac','MacStrict']:
  cols=np.flatnonzero(np.isin(g,strictgenes));take=np.ones(len(m),bool) if variant=='AllMac' else np.asarray((x[:,cols]>0).sum(1)).ravel()>=2
  for sam,mm in m.loc[take].groupby('sample'):
   ids=mm.index.to_numpy();lib=mm.library.sum();coverage.append(dict(dataset=ds,sample=sam,source=variant,cells=len(mm),library=lib,group=mm.group.iloc[0],day=int(mm.day.iloc[0])))
   for gene in needed:
    inds=np.flatnonzero(g==gene);v=np.asarray(x[ids][:,inds].sum(1)).ravel() if len(inds) else np.full(len(ids),np.nan)
    rows.append(dict(dataset=ds,sample=sam,group=mm.group.iloc[0],day=int(mm.day.iloc[0]),source=variant,gene=gene,cells=len(mm),UMI=v.sum(),CP10k=v.sum()/lib*1e4,fraction=float((v>0).mean()) if len(inds) else np.nan,feature_present=bool(len(inds))))
write(pd.DataFrame(rows),'sender_expression.tsv');write(pd.DataFrame(coverage),'sender_coverage.tsv')
sender=pd.DataFrame(rows);assert (pd.DataFrame(coverage).groupby('source').size()==10).all()
# Receiver readouts and backgrounds use complete measured genes, sample fractions.
rec=[];bgs=[]
for scope in ['Core','Naive','ActCD4','Treg']:
 cs=[];fs=[];met=[]
 for ds in ['GSE314943','GSE183682']:
  stem=ds+'__'+scope;c=pd.read_csv(B/'inputs'/(stem+'.counts.tsv.gz'),sep='\t',index_col=0);f=pd.read_csv(B/'inputs'/(stem+'.fraction.tsv.gz'),sep='\t',index_col=0);m=pd.read_csv(B/'inputs'/(stem+'.samples.tsv'),sep='\t')
  for r in m.itertuples():
   for gene in needed:
    rec.append(dict(dataset=ds,sample=r.sample,group=r.group,day=r.day,scope=scope,gene=gene,CP10k=float(c.at[gene,r.sample]/r.library*1e4) if gene in c.index else np.nan,fraction=float(f.at[gene,r.sample]) if gene in f.index else np.nan,cells=r.cells))
  cs.append(c);fs.append(f);met.append(m)
 common=sorted(set(cs[0].index)&set(cs[1].index));frac=pd.concat([f.loc[common] for f in fs],axis=1);eligible=(frac>=.1).sum(1)>=3
 for f in fs:eligible&=(f.loc[common]>=.1).sum(1)>=1
 bgs.extend(dict(scope=scope,gene=g) for g in frac.index[eligible])
write(pd.DataFrame(rec),'receiver_expression.tsv');write(pd.DataFrame(bgs),'joint_receiver_background.tsv')
rec=pd.DataFrame(rec);m=rec[rec.scope=='Core'][['dataset','sample','group','day']].drop_duplicates().sort_values(['dataset','sample']).reset_index(drop=True)
sampleids=list(zip(m.dataset,m['sample']));lrrows=[];lrsamples=[]
for source in ['AllMac','MacStrict']:
 s=sender[sender.source==source];scp=s.pivot(index='gene',columns=['dataset','sample'],values='CP10k').reindex(columns=sampleids);sfr=s.pivot(index='gene',columns=['dataset','sample'],values='fraction').reindex(columns=sampleids)
 for scope in ['Core','Naive','ActCD4','Treg']:
  r=rec[rec.scope==scope];rcp=r.pivot(index='gene',columns=['dataset','sample'],values='CP10k').reindex(columns=sampleids);rfr=r.pivot(index='gene',columns=['dataset','sample'],values='fraction').reindex(columns=sampleids)
  z=m.copy()
  for rr in lr.itertuples():
   lp=rr.ligand_parts.split(';');rp=rr.receptor_parts.split(';')
   measured=all(g in scp.index and scp.loc[g].notna().all() for g in lp) and all(g in rcp.index and rcp.loc[g].notna().all() for g in rp)
   if measured:
    sf=sfr.loc[lp].min(axis=0).to_numpy();rf=rfr.loc[rp].min(axis=0).to_numpy();sd=sf>=.1;rd=rf>=.1
    elig=sd.sum()>=3 and rd.sum()>=3 and all(sd[m.dataset.eq(ds)].any() and rd[m.dataset.eq(ds)].any() for ds in m.dataset.unique())
    lval=np.prod(scp.loc[lp].to_numpy(),axis=0)**(1/len(lp));rval=np.prod(rcp.loc[rp].to_numpy(),axis=0)**(1/len(rp));score=np.log1p(np.sqrt(lval*rval))
   else:elig=False;sd=rd=np.zeros(10,bool);lval=rval=score=np.full(10,np.nan)
   lrrows.append(dict(source=source,scope=scope,interaction=rr.interaction_name,ligand_parts=rr.ligand_parts,receptor_parts=rr.receptor_parts,annotation=rr.annotation,pathway=rr.pathway_name,measured=measured,eligible=bool(elig),sender_detected_n=int(sd.sum()),receiver_detected_n=int(rd.sum())))
   if elig:
    col='M_'+rr.interaction_name;z=m.assign(**{col:score})
    for i,(ds,sam) in enumerate(sampleids):lrsamples.append(dict(source=source,scope=scope,interaction=rr.interaction_name,dataset=ds,sample=sam,group=m.group.iloc[i],day=m.day.iloc[i],ligand_CP10k=lval[i],receptor_CP10k=rval[i],score=score[i]))
    S.group_effect(z,col,'Joint',source,scope)
write(pd.DataFrame(lrrows),'joint_LR_eligibility.tsv');write(pd.DataFrame(lrsamples),'joint_LR_sample_scores.tsv')
for row in S.effectrows:row['family']=row['family'].replace('programme_disease_effect','LR_disease_effect')
S.finalize(S.effectrows,'joint_LR_disease_effects.tsv')
# Whole-lung competing source expression is descriptive and separately labeled.
whole=[]
for ds,aa in [('GSE183682',a),('GSE314943',ad.read_h5ad(P/'02_analysis/GSE314943/stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad'))]:
 xx=sparse.csr_matrix(aa.layers['counts']);gg=np.array([canon(g) for g in (aa.var.symbol if 'symbol' in aa.var else aa.var_names)]);lib=np.asarray(xx.sum(1)).ravel()
 for (sam,broad),mm in aa.obs.groupby(['sample','Broad'],observed=True):
  ids=aa.obs_names.get_indexer(mm.index);zz=xx[ids];tot=lib[ids].sum()
  for gene in sorted(set(';'.join(lr.ligand_parts).split(';'))|set(cand.gene)):
   inds=np.flatnonzero(gg==gene)
   if not len(inds):continue
   vv=np.asarray(zz[:,inds].sum(1)).ravel();whole.append(dict(dataset=ds,sample=sam,Broad=broad,cells=len(ids),gene=gene,CP10k=vv.sum()/tot*1e4,fraction=float((vv>0).mean())))
write(pd.DataFrame(whole),'whole_lung_source_expression.tsv.gz')
(O/'logs/communication_export_complete.json').write_text(json.dumps({'status':'success','LR_eligible':int(pd.DataFrame(lrrows).eligible.sum()),'source_sample_rows':len(coverage)},indent=2),encoding='utf-8')
print('COMMUNICATION_COMPLETE',len(S.effectrows))
