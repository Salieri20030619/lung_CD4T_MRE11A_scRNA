from pathlib import Path
import pandas as pd,json,hashlib,shutil
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');B=P/'11_reports';S=B/'independent_candidates_revision_20261005_v1/results/source_evidence';O=B/'phase_delivery_IKZF2_STAT1_20261005_v1'
for s in ['results/source_evidence','provenance','code','deliverables','figures']:(O/s).mkdir(parents=True,exist_ok=True)
manifest=[];trace=[]
for tag in ['human_TF_disease','old_RNA_DE','old_TF_group','old_TF_program','IL7R_IKZF2_excess','old_M6_TF','mouse_TF_identity']:
 f=next(S.glob(tag+'.tsv*'));out=O/'results/source_evidence'/f.name;shutil.copy2(f,out);d=pd.read_csv(f,sep='\t',low_memory=False)
 manifest.append(dict(tag=tag,path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),rows=len(d)))
 pc=next((c for c in ['P','PValue','p','p_nominal'] if c in d),None)
 if pc:
  z=d.copy();z['source_tag']=tag;z['original_P_column']=pc;z['P_nominal_trace']=z[pc];z['nominal_P_lt_005']=z[pc].lt(.05);z['original_FDR_available']=z['FDR'].notna() if 'FDR' in z else False
  trace.append(z)
pd.concat(trace,ignore_index=True).to_csv(O/'results/TF_existing_models_full_trace_SERVER_ONLY.tsv.gz',sep='\t',index=False)
(O/'provenance/source_manifest.json').write_text(json.dumps(manifest,indent=2))
f=B/'aligned_evidence_20261002_v1/provenance/original_S2_analyze.py';ls=f.read_text().splitlines()
for i,l in enumerate(ls):
 if any(w in l for w in ['FDR','bh(','def assoc','SourceDep']):print(i+1,'\n'.join(ls[max(0,i-2):min(len(ls),i+3)]))
print('ARCHIVED',str(O),len(manifest),'source tables; no model rerun')
