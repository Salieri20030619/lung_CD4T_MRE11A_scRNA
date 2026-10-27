from pathlib import Path
import os,sys,json,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd,h5py
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/independent_candidates_revision_20261005_v1';G=P/'11_reports/GSE296682_patient_localization_20261005_v1';R=O/'results'
spec='Independent IL18 localization addendum. Frozen before extraction: IL18, IL18R1, IL18RAP; the same two official silicosis patients and all author subtypes plus AuthorMac_union. Stored count >0 detection; receptor same-cell joint detection. Existing TCR2/CD4_detected/MatureMac3 context only. No relabeling, new clustering, threshold optimization, DE, P or FDR. Cache-only official H5 read: fail on cache miss, no network. This fills an omitted candidate, not a candidate-to-candidate pathway.\n'
(O/'provenance/IL18_FINITE_SPEC.txt').write_text(spec)
sys.path.insert(0,str(G/'code'))
from g296_remote_h5 import OfficialRangeFile
class CacheOnly(OfficialRangeFile):
 def getblock(self,i):
  start=i*self.block;end=min(self.size-1,start+self.block-1);p=self.cache/f'{start}_{end}.bin'
  assert p.exists(),f'OFFLINE_CACHE_MISS: {p}'
  b=p.read_bytes();assert len(b)==end-start+1;return b
def dec(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
m=pd.read_csv(G/'results/patient_cell_marker_context_SERVER_ONLY.tsv.gz',sep='\t',index_col=0);genes=['IL18','IL18R1','IL18RAP'];f=CacheOnly()
with h5py.File(f,'r') as h:
 def index(group):
  key=group.attrs.get('_index','_index');key=key.decode() if isinstance(key,bytes) else key;return dec(group[key][:])
 names=index(h['var']);cells=index(h['obs']);lookup={x:i for i,x in enumerate(cells)};rows=[lookup[x] for x in m.index];cols=[int(np.flatnonzero(names==g)[0]) for g in genes];assert all((names==g).sum()==1 for g in genes);x=h['X'];enc=x.attrs['encoding-type'];enc=enc.decode() if isinstance(enc,bytes) else enc;assert enc=='csr_matrix';ptr=x['indptr'][:];X=np.zeros((len(m),3),np.float32)
 for k,r in enumerate(rows):
  lo,hi=int(ptr[r]),int(ptr[r+1]);ii=x['indices'][lo:hi];v=x['data'][lo:hi]
  for j,c in enumerate(cols):X[k,j]=v[ii==c].sum()
np.savez_compressed(R/'IL18_official_expression_SERVER_ONLY.npz',cells=np.array(m.index,dtype=str),genes=np.array(genes),values=X)
rows=[];context=[];scopes={s:m.Cell_Subtype.eq(s).to_numpy() for s in sorted(m.Cell_Subtype.unique())};scopes['AuthorMac_union']=m.Cell_Type.eq('Macrophage').to_numpy()
for sample in sorted(m.Sample.unique()):
 for scope,mask in scopes.items():
  sel=mask&m.Sample.eq(sample).to_numpy();N=int(sel.sum())
  for j,g in enumerate(genes):
   pos=sel&(X[:,j]>0);n=int(pos.sum());rows.append(dict(sample=sample,scope=scope,gene=g,N=N,n_positive=n,pct=n/N if N else np.nan,mean_stored=float(X[sel,j].mean()) if N else np.nan,status='NO_CELLS' if N==0 else 'ZERO_DETECTED' if n==0 else 'DETECTED'))
  joint=sel&(X[:,1]>0)&(X[:,2]>0);context.append(dict(sample=sample,scope=scope,N=N,receptor_joint=int(joint.sum()),receptor_joint_pct=float(joint.sum()/N) if N else np.nan,joint_TCR2=int((joint&m.TCR2.to_numpy()).sum()),joint_CD4=int((joint&m.CD4_detected.to_numpy()).sum()),joint_MatureMac3=int((joint&m.MatureMac3.to_numpy()).sum()),IL18_with_MatureMac3=int((sel&(X[:,0]>0)&m.MatureMac3.to_numpy()).sum())))
q=pd.DataFrame(rows);q.to_csv(R/'IL18_localization_all_author_scopes.tsv',sep='\t',index=False);c=pd.DataFrame(context);c.to_csv(R/'IL18_receptor_joint_context.tsv',sep='\t',index=False)
# Reused network input audit: positive importance weights do not supply a calibrated signed intervention model.
n=pd.read_csv(R/'source_evidence/network_edges.tsv',sep='\t');a=n.groupby(['species','filter','seed','tf']).agg(edges=('target','size'),targets=('target','nunique'),min_weight=('weight','min'),max_weight=('weight','max'),nonpositive=('weight',lambda s:int((s<=0).sum())),retained=('retained_DropOverlap','sum')).reset_index();a['virtual_KO_status']='FROZEN_REGULON_IMPORTANCE_AVAILABLE;SIGNED_INTERVENTION_COEFFICIENT_AND_VALIDATION_NOT_ESTABLISHED';a.to_csv(R/'TF_virtual_KO_input_audit.tsv',sep='\t',index=False)
info=dict(status='COMPLETE',genes=genes,cells=len(m),patients=2,new_network_requests=f.requests,whole_file_checksum_verified=False,source='existing_official_Zenodo_range_cache',P_FDR='NOT_PERFORMED_LOCALIZATION_ONLY',spec_sha256=hashlib.sha256(spec.encode()).hexdigest());(O/'provenance/IL18_addendum_complete.json').write_text(json.dumps(info,indent=2));focus=['Alveolar macrophage','Monocyte-like macrophage','Interstitial-like macrophage','CD4 effector T cell','Regulatory T cell'];print(q[q.scope.isin(focus)].to_string(index=False));print(c[c.scope.isin(focus)].to_string(index=False));print(a.to_string(index=False));print(json.dumps(info),flush=True)
