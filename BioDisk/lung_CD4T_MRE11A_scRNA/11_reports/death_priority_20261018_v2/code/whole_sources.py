import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
from pathlib import Path
import anndata as ad,pandas as pd,numpy as np,json
from scipy import sparse
O=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v2');P=O.parent.parent;B=P/'11_reports/three_cohort_four_axis_20261017_v16'
lr=pd.read_csv(B/'resources/CellChat_protein_LR_expanded.tsv',sep='\t').fillna('');cand=pd.read_csv(O/'tables/frozen_candidate_pairs.tsv',sep='\t')
genes=sorted((set(';'.join(lr.ligand_parts).split(';'))|set(cand.gene))-{''})
alias={'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}
rows=[]
for ds,path in [('GSE183682',O/'objects/GSE183682_whole_clustered.h5ad'),('GSE314943',P/'02_analysis/GSE314943/stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad')]:
 a=ad.read_h5ad(path)
 if ds=='GSE183682':a.obs['Broad']=a.obs.cluster.astype(str).map(json.loads((O/'inputs/GSE183682_annotation.json').read_text()))
 x=sparse.csr_matrix(a.layers['counts']);g=np.array([alias.get(str(v).upper(),str(v).upper()) for v in (a.var.symbol if 'symbol' in a.var else a.var_names)]);lib=np.asarray(x.sum(1)).ravel()
 # Collapse mapped features once, slice the selected small matrix before group work.
 lookup={v:i for i,v in enumerate(genes)};ii=np.flatnonzero(np.isin(g,genes));proj=sparse.csr_matrix((np.ones(len(ii)),(ii,[lookup[g[i]] for i in ii])),shape=(len(g),len(genes)));selected=(x@proj).tocsr();present=set(g)
 for (sam,broad),mm in a.obs.groupby(['sample','Broad'],observed=True):
  ids=a.obs_names.get_indexer(mm.index);z=selected[ids];cp=np.asarray(z.sum(0)).ravel()/lib[ids].sum()*1e4;fraction=np.asarray((z>0).sum(0)).ravel()/len(ids)
  rows.extend(dict(dataset=ds,sample=sam,Broad=broad,cells=len(ids),gene=gene,CP10k=cp[j],fraction=fraction[j]) for j,gene in enumerate(genes) if gene in present)
 print('SOURCE',ds,len(rows),flush=True)
pd.DataFrame(rows).to_csv(O/'tables/whole_lung_source_expression.tsv.gz',sep='\t',index=False)
eligible=pd.read_csv(O/'tables/joint_LR_eligibility.tsv',sep='\t');cover=pd.read_csv(O/'tables/sender_coverage.tsv',sep='\t')
(O/'logs/communication_export_complete.json').write_text(json.dumps({'status':'success','LR_eligible':int(eligible.eligible.sum()),'source_sample_rows':len(cover),'source_export':'optimized equivalent sparse projection; original slow source loop terminated after LR outputs completed'},indent=2),encoding='utf-8')
print('COMMUNICATION_COMPLETE')
