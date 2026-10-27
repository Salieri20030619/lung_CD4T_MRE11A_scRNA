from pathlib import Path
import json,re,hashlib
import numpy as np,pandas as pd
from scipy import sparse,io
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018')
B=P/'11_reports/three_cohort_four_axis_20261017_v16/inputs'
V=P/'11_reports/death_priority_20261018_v2'
xs=[];gs=[];ms=[]
for ds in ['GSE314943','GSE183682']:
 x=sparse.load_npz(B/(ds+'__Core.cellcounts.npz')).tocsr();g=pd.read_csv(B/(ds+'__Core.genes.txt'),header=None)[0].to_numpy();m=pd.read_csv(B/(ds+'__Core.cells.tsv.gz'),sep='\t');m['dataset']=ds;m['cell']=ds+'|'+m.cell
 assert x.shape==(len(m),len(g)) and np.all(x.data>=0) and np.all(x.data==np.floor(x.data))
 xs.append(x);gs.append(g);ms.append(m)
common=sorted(set(gs[0])&set(gs[1]));x=sparse.vstack([a[:,pd.Index(g).get_indexer(common)] for a,g in zip(xs,gs)],format='csr');m=pd.concat(ms,ignore_index=True);g=np.array(common)
modules=pd.read_csv(V/'tables/module_members.tsv',sep='\t');modules.to_csv(O/'module_members.tsv',sep='\t',index=False)
print('MODULE COLUMNS',modules.columns.tolist(),flush=True)
program=set(modules.gene)
excluded=np.array([bool(re.match(r'^(RP[LS][0-9]|RPLP[0-9]|RPSA$|MRP[LS][0-9]|MT-)',s)) for s in g])
rows=[];jobs=[];geneaudit=[]
for arm in ['All','BG1']:
 for sub,targets in [('Naive',['HIF1A','ARL6IP1']),('Treg',['IL2RA'])]:
  for state in ['Saline','Silica']:
   mask=m.subtype.eq(sub)&m.group.eq(state)
   if arm=='BG1':mask &= m.background_fraction.le(.01)
   inds=np.flatnonzero(mask);full=x[inds];det=np.asarray((full>0).mean(axis=0)).ravel();lib=np.asarray(full.sum(axis=1)).ravel();norm=full.multiply(1e4/lib[:,None]).tocsr();norm.data=np.log1p(norm.data);mean=np.asarray(norm.mean(axis=0)).ravel();var=np.asarray(norm.power(2).mean(axis=0)).ravel()-mean**2
   eligible=(det>=.05)&(var>1e-8)&~excluded
   # Select variable genes within mean-expression bins, without disease comparison.
   score=np.full(len(g),-np.inf);ii=np.flatnonzero(eligible);bins=pd.qcut(mean[ii],20,labels=False,duplicates='drop')
   for bn in np.unique(bins):
    jj=ii[np.asarray(bins)==bn];v=np.log(var[jj]+1e-8);score[jj]=(v-v.mean())/(v.std()+1e-8)
   chosen=set(g[np.argsort(-score)[:min(3000,eligible.sum())]])|((program|set(targets))&set(g[eligible]));chosen=sorted(chosen);jj=pd.Index(g).get_indexer(chosen)
   for gene in targets:
    k=pd.Index(g).get_loc(gene);geneaudit.append(dict(arm=arm,subtype=sub,state=state,gene=gene,cells=len(inds),detection=det[k],eligible=bool(eligible[k]),retained=gene in chosen))
   for seed in [17,29,43]:
    rng=np.random.default_rng(seed);selected=[]
    for sam,mm in m.loc[inds].groupby('sample'):
     take=rng.choice(mm.index.to_numpy(),min(150,len(mm)),replace=False);selected.extend(take.tolist())
    selected=sorted(selected);out=O/'inputs'/f'{arm}_{sub}_{state}_s{seed}';out.mkdir(parents=True,exist_ok=True)
    xx=x[selected][:,jj].T.tocoo();io.mmwrite(out/'counts.mtx',xx,field='integer');(out/'genes.txt').write_text('\n'.join(chosen)+'\n',encoding='utf-8');m.loc[selected].to_csv(out/'cells.tsv',sep='\t',index=False)
    avail=[t for t in targets if t in chosen]
    job=dict(id=out.name,path=str(out),arm=arm,subtype=sub,state=state,seed=seed,targets=avail,genes=len(chosen),cells=len(selected),donors=int(m.loc[selected,'sample'].nunique()))
    jobs.append(job)
    for sam,mm in m.loc[selected].groupby('sample'):rows.append(dict(job=out.name,sample=sam,dataset=mm.dataset.iloc[0],day=int(mm.day.iloc[0]),cells=len(mm)))
pd.DataFrame(rows).to_csv(O/'input_donors.tsv',sep='\t',index=False);pd.DataFrame(geneaudit).to_csv(O/'target_coverage.tsv',sep='\t',index=False)
(O/'jobs.json').write_text(json.dumps(jobs,indent=2),encoding='utf-8')
(O/'protocol.json').write_text(json.dumps(dict(method='scTenifoldKnk 2.0.0 classic manifold, independent single gene knockouts',seeds=[17,29,43],arms=['All','BG1 <=0.01 ambient fraction'],states='Saline and Silica fit separately; WT means target intact, not necessarily disease-free',genes='top 3000 mean-binned logCP10k variance plus eligible prespecified program/target genes, detection >=5%; ribosomal and mitochondrial genes excluded',sampling='cap 150 cells per donor without replacement per seed; max500 cells/network,10 networks',statistics='model P and BH within each knockout/run, not biological treatment P; gene ranks and pathway overlap across seeds, no pooling cells as biological replicates',direction='secondary diffusion prediction with library-size regression, distinct from unsigned manifold distance',limitations='controls Treg only125 cells; cohort/day imbalance remains; no claim of exact biological KO or cell death'),indent=2),encoding='utf-8')
print(pd.DataFrame(jobs).to_string(index=False),flush=True)
print(pd.DataFrame(geneaudit).to_string(index=False),flush=True)