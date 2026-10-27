from pathlib import Path
import sys,json,hashlib,shutil
import numpy as np,pandas as pd,h5py
O=Path(sys.argv[1]);P=O.parents[1];J=P/'11_reports/GSE296682_patient_localization_20261005_v1';(O/'public_cache').mkdir(exist_ok=True);cache=O/'public_cache/h5ad_range_blocks';cache.mkdir(exist_ok=True)
for f in (J/'public_cache/h5ad_range_blocks').glob('*.bin'):
 t=cache/f.name
 if not t.exists():t.hardlink_to(f)
f=J/'public_cache/BAL_FINAL.h5ad.range0_65535'
if f.exists() and not (O/'public_cache'/f.name).exists():(O/'public_cache'/f.name).hardlink_to(f)
ns={'__name__':'range_reuse'};exec((J/'code/g296_remote_h5.py').read_text().replace('96*1024**2','128*1024**2'),ns);ns['O']=O;stream=ns['OfficialRangeFile']();spec=json.loads((O/'provenance/specification.json').read_text());G=spec['genes'];pairs=spec['complexes'];m=pd.read_csv(J/'results/official_silicosis_metadata_SERVER_ONLY.tsv.gz',sep='\t',index_col=0)
def dec(a):return np.array([z.decode() if isinstance(z,bytes) else str(z) for z in a])
with h5py.File(stream,'r') as h:
 cells=dec(h['obs'][h['obs'].attrs['_index']][:]);gn=dec(h['var'][h['var'].attrs['_index']][:]);rowmap={g:i for i,g in enumerate(cells)};assert set(m.index)<=set(rowmap);rr=np.array([rowmap[g] for g in m.index]);look=np.full(len(cells),-1);look[rr]=np.arange(len(rr));x=h['X'];encoding=x.attrs['encoding-type'];ptr=x['indptr'][:];data=np.full((len(rr),len(G)),np.nan);mapping=[]
 ixall=[]
 for j,g in enumerate(G):
  ix=np.flatnonzero(gn==g);mapping.append(dict(gene=g,measured=len(ix)==1));ixall.append(int(ix[0]) if len(ix)==1 else -1)
 if encoding in ['csc_matrix',b'csc_matrix']:
  for j,col in enumerate(ixall):
   if col<0:continue
   lo,hi=int(ptr[col]),int(ptr[col+1]);ii=x['indices'][lo:hi];vv=x['data'][lo:hi];take=look[ii]>=0;data[:,j]=0;data[look[ii[take]],j]=vv[take]
 else:
  assert encoding in ['csr_matrix',b'csr_matrix'];cols=np.full(len(gn),-1);valid=np.array(ixall)>=0;cols[np.array(ixall)[valid]]=np.flatnonzero(valid);data[:,valid]=0
  for k,row in enumerate(rr):
   lo,hi=int(ptr[row]),int(ptr[row+1]);ii=x['indices'][lo:hi];vv=x['data'][lo:hi];take=cols[ii]>=0;data[k,cols[ii[take]]]=vv[take]
   if k%1000==0:print('ROWS',k,len(rr),flush=True)
rows=[];joint=[];scope_masks={}
for label in sorted(m.Cell_Subtype.dropna().unique()):scope_masks[label]=m.Cell_Subtype.eq(label).to_numpy()
for label,mask0 in scope_masks.items():
 for sample in sorted(m.Sample.unique()):
  mask=mask0&m.Sample.eq(sample).to_numpy();n=int(mask.sum())
  for j,g in enumerate(G):
   vals=data[mask,j];ok=n>0 and np.isfinite(vals).all();rows.append(dict(sample=sample,scope=label,gene=g,cells=n,measured=bool(np.isfinite(data[:,j]).all()),positive_cells=int((vals>0).sum()) if ok else np.nan,pct=float((vals>0).mean()*100) if ok else np.nan,mean_stored=float(vals.mean()) if ok else np.nan))
  for name,gs in pairs.items():
   ix=[G.index(g) for g in gs];z=data[mask][:,ix];ok=n>0 and np.isfinite(z).all();pos=(z>0).all(axis=1);joint.append(dict(sample=sample,scope=label,complex=name,components=';'.join(gs),cells=n,positive_cells=int(pos.sum()) if ok else np.nan,pct=float(pos.mean()*100) if ok else np.nan))
pd.DataFrame(rows).to_csv(O/'tables/patient2_localization.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(joint).to_csv(O/'tables/patient2_complex_codetection.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(mapping).to_csv(O/'tables/patient2_gene_coverage.tsv',sep='\t',index=False)
np.savez_compressed(O/'tables/patient2_selected_expression.npz',values=data,cells=np.array(m.index,dtype=str),genes=np.array(G,dtype=str))
(O/'provenance/patient2_provenance.json').write_text(json.dumps({'source':stream.url,'selected_patients':sorted(m.Sample.unique()),'genes':G,'new_range_requests':stream.requests,'source_full_md5':'86e4165b84f965e807e4cf361d3d6bb6','whole_file_checksum_verified':False,'interpretation':'patient-only localization, stored processed expression, no disease effects or n=2 correlations','metadata_sha256':hashlib.sha256((J/'results/official_silicosis_metadata_SERVER_ONLY.tsv.gz').read_bytes()).hexdigest()},indent=2),encoding='utf-8');print('DONE',len(m),len(rows),stream.requests)