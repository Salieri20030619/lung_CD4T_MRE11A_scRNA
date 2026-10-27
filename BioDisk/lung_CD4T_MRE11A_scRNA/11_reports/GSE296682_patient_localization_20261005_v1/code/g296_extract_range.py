from pathlib import Path
import os,json,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd,h5py
from scipy import sparse
from g296_remote_h5 import OfficialRangeFile,O
R=O/'results';m=pd.read_csv(R/'official_silicosis_metadata_SERVER_ONLY.tsv.gz',sep='\t',index_col=0)
fixed='IL15 TNF CLSTN1 ADAM15 PLAU DSC2 CXCL2 ICAM2 VTN VEGFB GSTP1 RTN4 IL15RA IL2RB IL2RG TNFRSF1A TNFRSF1B IKZF2 STAT1'.split()
panels={'T':'CD3D CD3E TRAC CD4 IL7R LTB CCR7 SELL TCF7','CD8_NK':'CD8A CD8B NKG7 GNLY KLRD1 PRF1','Treg':'FOXP3 IL2RA CTLA4 TIGIT','Myeloid':'LYZ LST1 TYROBP AIF1 FCER1G','MatureMac':'C1QA C1QB C1QC APOE APOC1 MRC1 MSR1 CD68 CD163 CSF1R MERTK TREM2 GPNMB LGMN SPP1 FOLR2 LYVE1','AM':'FABP4 PPARG MARCO','Mono':'FCN1 VCAN S100A8 S100A9','Competitor':'CD1C FCER1A CLEC10A EPCAM KRT8 KRT18 SCGB1A1 PPBP PF4'}
genes=list(dict.fromkeys(fixed+sum([v.split() for v in panels.values()],[])))
(O/'provenance/fixed_genes_and_review_panels.json').write_text(json.dumps(dict(candidates=fixed,panels=panels),indent=2))
def dec(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
f=OfficialRangeFile()
with h5py.File(f,'r') as h:
 var=h['var'];vindex=var.attrs.get('_index','_index');vindex=vindex.decode() if isinstance(vindex,bytes) else vindex;names=dec(var[vindex][:]);obs=h['obs'];oid=obs.attrs.get('_index','_index');oid=oid.decode() if isinstance(oid,bytes) else oid;cells=dec(obs[oid][:]);assert len(set(cells))==len(cells);assert len(set(names))==len(names)
 rowmap={x:i for i,x in enumerate(cells)};assert set(m.index)<=set(rowmap);rows=np.array([rowmap[x] for x in m.index]);mapping=[];colids=[];measured=[]
 for g in genes:
  ix=np.flatnonzero(names==g)
  if len(ix)==0:ix=np.flatnonzero(np.char.upper(names)==g.upper())
  ok=len(ix)==1;mapping.append(dict(gene=g,status='EXACT' if ok and names[ix[0]]==g else 'UNIQUE_CASE' if ok else 'ABSENT_OR_AMBIGUOUS',feature=names[ix[0]] if ok else '',column=int(ix[0]) if ok else -1,stable_id_available=False))
  if ok:colids.append(int(ix[0]));measured.append(g)
 pd.DataFrame(mapping).to_csv(R/'fixed_feature_mapping.tsv',sep='\t',index=False)
 x=h['X'];encoding=x.attrs.get('encoding-type','');encoding=encoding.decode() if isinstance(encoding,bytes) else encoding;shape=list(map(int,x.attrs['shape']));print('EXTRACT_LAYOUT',encoding,shape,'selected_cells',len(rows),'genes',len(measured),'var_columns',list(var.keys()),'layers',list(h['layers'].keys()),'raw_present','raw' in h,flush=True)
 assert shape==[len(cells),len(names)];data=np.zeros((len(rows),len(colids)),dtype=np.float32);ptr=x['indptr'][:]
 if encoding=='csc_matrix':
  rowlookup=np.full(len(cells),-1,dtype=np.int32);rowlookup[rows]=np.arange(len(rows))
  for done,jj in enumerate(np.argsort(colids)):
   col=colids[jj];lo,hi=int(ptr[col]),int(ptr[col+1]);ix=x['indices'][lo:hi];val=x['data'][lo:hi];take=rowlookup[ix]>=0;data[rowlookup[ix[take]],jj]=val[take]
   if done%10==0:print('EXTRACT_GENE_PROGRESS',done+1,len(colids),measured[jj],flush=True)
 elif encoding=='csr_matrix':
  colslookup=np.full(len(names),-1,dtype=np.int32);colslookup[colids]=np.arange(len(colids))
  for k,row in enumerate(rows):
   lo,hi=int(ptr[row]),int(ptr[row+1]);ix=x['indices'][lo:hi];val=x['data'][lo:hi];take=colslookup[ix]>=0;data[k,colslookup[ix[take]]]=val[take]
   if k%500==0:print('EXTRACT_ROW_PROGRESS',k,len(rows),flush=True)
 else:raise RuntimeError('Unexpected matrix storage, stop for review')
 info=dict(source='OFFICIAL_ZENODO_V2_RANGE_READ',url=f.url,official_full_bytes=f.size,official_full_md5='86e4165b84f965e807e4cf361d3d6bb6',whole_file_checksum_verified=False,range_blocks_hash_log=str(O/'provenance/range_blocks.jsonl'),shape=shape,X_encoding=encoding,layers=list(h['layers'].keys()),raw_present='raw' in h,var_columns=list(var.keys()),selected_cells=len(rows),selected_genes=len(measured),all_values_nonnegative=bool((data>=0).all()),positive_values_integer_fraction=float(np.mean(np.isclose(data[data>0],np.round(data[data>0]),rtol=0,atol=1e-6))),max=float(data.max()))
np.savez_compressed(R/'official_fixed_expression_SERVER_ONLY.npz',values=data,cells=np.array(m.index,dtype=str),genes=np.array(measured));(O/'provenance/official_fixed_expression_provenance.json').write_text(json.dumps(info,indent=2));print('OFFICIAL_EXTRACTION_COMPLETE',json.dumps(info),flush=True)
# Compare only the targeted cells/genes to cached GEO, never assume all matrices equivalent.
with h5py.File('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/00_raw/GSE296682/GSE296682_BAL_FINAL.h5ad','r') as old:
 on=dec(old['var/_index'][:]);oc=dec(old['obs/_index'][:]);ori={x:i for i,x in enumerate(oc)};orows=np.array([ori[x] for x in m.index]);ox=old['X'];op=ox['indptr'][:];lookup=np.full(len(oc),-1,dtype=np.int32);lookup[orows]=np.arange(len(orows));comp=[]
 for j,g in enumerate(measured):
  ix=np.flatnonzero(on==g)
  if len(ix)!=1:comp.append(dict(gene=g,status='OLD_FEATURE_UNRESOLVED'));continue
  col=ix[0];lo,hi=int(op[col]),int(op[col+1]);ii=ox['indices'][lo:hi];v=ox['data'][lo:hi];q=np.zeros(len(orows));take=lookup[ii]>=0;q[lookup[ii[take]]]=v[take];comp.append(dict(gene=g,status='COMPARED',cells=len(q),exact=int((q==data[:,j]).sum()),different=int((q!=data[:,j]).sum()),max_abs_diff=float(np.max(abs(q-data[:,j])))))
pd.DataFrame(comp).to_csv(R/'fixed_expression_GEO_Zenodo_comparison.tsv',sep='\t',index=False);print('CACHE_COMPARISON',pd.DataFrame(comp).fillna(0).different.sum(),flush=True)
