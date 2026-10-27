import os,json,hashlib,shutil,tarfile,sys
from pathlib import Path
import numpy as np,pandas as pd
import loompy
from pyscenic.cli.utils import load_signatures
ROOT=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
OUT=ROOT/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1'
OUT.mkdir(parents=True,exist_ok=False)
for d in ['source_tables','existing_figures','provenance','data','code']:(OUT/d).mkdir()
S1=ROOT/'02_analysis/scenic1_20260929_103513_xRoEnr'
S2=ROOT/'02_analysis/scenic2_20261001_092655_ZKNWD6'
S3=ROOT/'02_analysis/scenic3_20261001_123452_Kvp9h5'
S4=ROOT/'02_analysis/scenic4_20261002_091545_iazXf1'
M={1:ROOT/'02_analysis/mechanism_step1_20260927_152218_s4V5eF',2:ROOT/'02_analysis/mechanism_step2_20260927_154919_EsN4gk',3:ROOT/'02_analysis/mechanism_step3_20260928_103057_lIRL38',4:ROOT/'02_analysis/mechanism_step4_20260928_150337_J58beC',5:ROOT/'02_analysis/mechanism_step5_20260928_153825_Q4sB1x',6:ROOT/'02_analysis/mechanism_step6_20260928_163015_VoRQqA',7:ROOT/'02_analysis/mechanism_step7_20260929_100321_6hm4dZ',8:ROOT/'02_analysis/mechanism_step8_20260929_101407_jWFgLl'}
sources={}
def record(p):
 p=Path(p)
 if str(p) not in sources:
  st=p.stat()
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
  sources[str(p)]=dict(source=str(p),bytes=st.st_size,mtime_ns=st.st_mtime_ns,sha256=h.hexdigest())
 return p
def rd(p,**kw):return pd.read_csv(record(p),sep='\t',**kw)
def save(df,name):df.to_csv(OUT/'data'/name,sep='\t',index=False)
def copy(p,name,dest='source_tables'):shutil.copy2(record(p),OUT/dest/name)
def canon(g):
 s=str(g).upper()
 return {'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'}.get(s,s)
for n,base in [('M'+str(i),p) for i,p in M.items()]+[('S'+str(i),p) for i,p in enumerate([S1,S2,S3,S4],1)]:
 if (base/'methods.txt').exists():copy(base/'methods.txt',n+'_methods.txt','provenance')
 if (base/'status.json').exists():copy(base/'status.json',n+'_status.json','provenance')
for name in ['TF_group_effects.tsv','TF_associations.tsv','TF_sample_activity.tsv.gz','regulon_exclusion_coverage.tsv']:
 d=rd(S2/'tables'/name)
 field='regulon' if 'regulon' in d else 'tf'
 keep=d[field].astype(str).str.split('(').str[0].str.upper().isin(['IKZF2','STAT1'])
 d=d[keep].copy()
 d['tf']=d[field].astype(str).str.split('(').str[0].str.upper()
 if name=='TF_associations.tsv':d=d[d.feature.str.startswith('score_')]
 save(d,'S2_'+name.replace('.gz',''))
print('FOCUSED_SCENIC_TABLES_EXPORTED',flush=True)
for n in ['mouse_TF_identity_adjustment.tsv','mouse_marker_sample.tsv','human_cluster_scorecard.tsv','human_cluster_markers.tsv','human_cluster_marker_expression.tsv','human_cluster_sample_composition.tsv','human_marker_coverage.tsv','human_sender_cluster_expression.tsv','human_sender_parent_sensitivity.tsv','human_sender_sample_covariates.tsv','CD4_receptor_detection.tsv','CD4_receptor_codetection.tsv']:
 copy(S3/'tables'/n,'S3_'+n)
for n in ['external_time_matched_TF_deltas.tsv','direction_replication.tsv','transfer_coverage.tsv','external_time_matched_phenotype_deltas.tsv','input_design.tsv']:
 copy(S4/'tables'/n,'S4_'+n)
for prefix,p,files in [
 ('M1',M[1],['module_definition.tsv','module_coverage.tsv','group_effects.tsv','sample_scores.tsv']),
 ('M2',M[2],['human_subject_design.tsv','human_source_models.tsv','human_sample_scores.tsv','human_module_coverage.tsv']),
 ('M3',M[3],['geneset_status.tsv','ligand_activities.tsv','supported_LR_edges.tsv','hypothesis_network.tsv']),
 ('M5',M[5],['focused_ligand_ranks.tsv','geneset_status.tsv']),
 ('M8',M[8],['annotation_counts.tsv','annotation_marker_summary.tsv','sample_mapping_status.tsv','sample_mapping_candidates.tsv','sample_mapping_edges.tsv','CD4_annotation_crosscheck.tsv'])]:
 for n in files:copy(p/'tables'/n,prefix+'_'+n)
e=rd(M[6]/'tables/TF_evidence_matrix.tsv')
save(e[e.tf.str.upper().isin(['IKZF2','STAT1'])],'M6_focused_TF_evidence.tsv')
de=[]
for arm in ['All','BG1']:
 for scope in ['Core','ActCD4','Naive','Treg']:
  p=M[3]/'tables'/f'DE_CD4_{scope}_{arm}.tsv'
  if p.exists():
   z=rd(p);z=z[z.gene.map(canon).isin(['IKZF2','STAT1'])].copy()
   z['species']='Mouse';z['scope']=scope;z['filter']=arm;z['source_file']=str(p);z['model']='edgeR_QL_group'
   de.append(z)
p=ROOT/'10_tables/CD4_final/pathway_v1/source_adjusted_full_DEG.tsv'
z=rd(p);z=z[z.gene.map(canon).isin(['IKZF2','STAT1'])].copy();z['species']='Human';z['scope']=z.analysis;z['filter']='All';z['source_file']=str(p);z['model']='edgeR_source_adjusted';de.append(z)
de=pd.concat(de,ignore_index=True);de['tf']=de.gene.map(canon);save(de,'existing_RNA_DE_focused.tsv')
for p in (ROOT/'10_tables/CD4_final').glob('*.tsv'):
 if any(x in p.name for x in ['composition','proportion']):copy(p,'Human_'+p.name)
copy(ROOT/'10_tables/broad_clustering_v1/GSE288887_preliminary_lineage_summary.tsv','Human_preliminary_lineage_summary.tsv')
mouse=ROOT/'02_analysis/GSE314943/stage5_20260905_222457_a2dApT'
for n in ['subtype_percent.tsv','subtype_counts.tsv','sample_exploratory_tests.tsv','PB_sample_metadata.tsv']:
 copy(mouse/'tables'/n,'Mouse_'+n)
copy(mouse/'tables/methods.txt','Mouse_CD4_methods.txt','provenance')
figs=[
(ROOT/'09_figures/broad_clustering_v1/GSE288887_preliminary_lineages.png','human_preliminary_atlas.png'),
(ROOT/'09_figures/broad_clustering_v1/GSE288887_lineage_marker_dotplot.png','human_lineage_markers.png'),
(ROOT/'09_figures/CD4_final/GSE288887_CD4_final_UMAP.png','human_CD4_UMAP.png'),
(ROOT/'09_figures/CD4_final/robustness/CD4_subject_composition_corrected.png','human_CD4_composition.png'),
(ROOT/'02_analysis/GSE314943/stage3_20260905_214228_VPtd0r/figures/whole_Broad.png','mouse_whole_atlas.png'),
(mouse/'figures/CD4_UMAP.png','mouse_CD4_UMAP.png'),(mouse/'figures/CD4_comp.png','mouse_CD4_composition.png'),
(M[3]/'figures/NicheNet_Core_Up.png','NicheNet_Core_Up.png'),(M[3]/'figures/NicheNet_Core_Down.png','NicheNet_Core_Down.png'),
(M[1]/'figures/module_scores.png','M1_module_scores.png'),(M[2]/'figures/human_modules.png','M2_human_modules.png'),
(S3/'figures/MacClust.png','S3_MacClust.png'),(S4/'figures/ExtTF.png','S4_ExtTF.png')]
for p,n in figs:copy(p,n,'existing_figures')
mods={canon(g) for g in rd(M[1]/'tables/module_definition.tsv').gene}
old={'IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11'}
rna=[];edges=[]
for species in ['Mouse','Human']:
 for arm in ['All','BG1']:
  base=S1/species/arm
  meta=rd(base/'cells.tsv')
  lp=record(base/'score.loom')
  with loompy.connect(str(lp),mode='r') as ds:
   genes=[canon(x) for x in ds.ra.Gene];cells=pd.Index(ds.ca.CellID.astype(str))
   assert cells.is_unique and set(cells)==set(meta.cell)
   idx=cells.get_indexer(meta.cell)
   for tf in ['IKZF2','STAT1']:
    present=tf in genes
    val=np.expm1(ds[genes.index(tf),:].reshape(-1)[idx]).astype(float) if present else np.full(len(meta),np.nan)
    for scope in ['Core']+sorted(meta.subtype.unique()):
     m=meta if scope=='Core' else meta[meta.subtype==scope]
     for sample,q in m.groupby('sample'):
      v=val[q.index.to_numpy()]
      rna.append(dict(species=species,filter=arm,scope=scope,sample=sample,group=q.group.iloc[0],series=q.series.iloc[0],tf=tf,cells=len(q),eligible20=len(q)>=20,present=present,mean_CP10k=float(np.mean(v)),pct=float(100*np.mean(v>0)) if present else np.nan,source=str(lp),transformation='expm1(frozen log1p(original-library CP10k)); sample arithmetic mean; no new statistical fit'))
  for seed in [17,29,43]:
   rp=record(base/f'seed{seed}/regulons.csv')
   for r in load_signatures(str(rp)):
    tf=canon(r.transcription_factor)
    if tf not in ['IKZF2','STAT1']:continue
    for g,w in r.gene2weight.items():
     target=canon(g)
     edges.append(dict(species=species,filter=arm,seed=seed,tf=tf,target=target,weight=float(w),is_self=target==tf,is_old_candidate=target in old,is_module_gene=target in mods,retained_DropOverlap=target not in old|mods|{tf},source=str(rp)))
  print('RNA_AND_NETWORK_EXPORTED='+species+' '+arm,flush=True)
save(pd.DataFrame(rna),'RNA_sample_summary.tsv');save(pd.DataFrame(edges),'TF_network_all_edges.tsv')
# Explicit audit of the user-supplied historical working mapping; no new annotation inference.
sc=rd(S3/'tables/human_cluster_scorecard.tsv');sc=sc[sc.resolution=='r4'].copy()
mix=rd(S3/'tables/human_cluster_sample_composition.tsv');mix=mix[mix.resolution=='r4'].copy();mix['cluster']=mix.cluster.astype(str);mix['parent']=mix.parent.astype(str)
mk=rd(S3/'tables/human_cluster_markers.tsv');mk['group']=mk['group'].astype(str)
working={0:'MacAM',2:'MacAM',11:'MacAM',13:'MacAM',3:'MonoMac',7:'MonoMac',1:'Mono-like',4:'Mono-like',5:'Mono-like',14:'Mono-like',6:'DC2-like',12:'DC1-like',8:'Epithelial-myeloid mixed',9:'T-myeloid mixed',10:'Vascular/epithelial/myeloid mixed'}
rows=[]
for _,s in sc.iterrows():
 c=str(s.cluster);v=mix[mix.cluster==c];total=int(v.cells.sum())
 bydonor=v[v['sample']!='Unknown'].groupby('sample').cells.sum()
 bysource=v.groupby('series').cells.sum()
 par=v.groupby('parent').cells.sum()
 top=mk[mk['group']==c].sort_values('scores',ascending=False).head(10)['names'].tolist()
 row=dict(cluster=c,working_label_from_prior_review=working.get(int(c),'Unresolved'),approval_status='Previously delivered working annotation; not per-label user-finalized',cells=total,mapped_subjects=len(bydonor),eligible20_subjects=int((bydonor>=20).sum()),unknown_cells=int(v.loc[v['sample']=='Unknown','cells'].sum()),dominant_mapped_subject=bydonor.idxmax() if len(bydonor) else '',dominant_subject_fraction=float(bydonor.max()/total) if len(bydonor) else np.nan,parent16_fraction=float(par.get('16',0)/total),parent_composition='; '.join(f'{k}:{int(n)}' for k,n in par.items()),source_composition='; '.join(f'{k}:{int(n)}' for k,n in bysource.items()),program_rank1=s.program_rank1,program_rank2=s.program_rank2,epi_fraction_mean=s.epi_fraction_mean,top10_markers=';'.join(top),review_action='Reconcile against original reviewed_annotation.tsv when supplied; no relabeling or model fit performed')
 rows.append(row)
review=pd.DataFrame(rows);save(review,'myeloid_working_annotation_audit.tsv')
print('WORKING_MAPPING_COUNTS='+repr(review.groupby('working_label_from_prior_review').cells.sum().to_dict()),flush=True)
manifest=pd.DataFrame(sources.values());manifest.to_csv(OUT/'provenance/source_manifest.tsv',sep='\t',index=False)
(OUT/'provenance/extraction_notes.txt').write_text('Read-only extraction from frozen outputs. Original files unchanged. New statistical models, thresholds, clustering, annotation decisions and GRNs were NOT run. RNA_sample_summary inverts the documented stored log1p CP10k and aggregates by existing sample/scope; existing_RNA_DE_focused retains original edgeR P/FDR and scope. Historical working annotation was supplied by parent conversation and awaits original review ZIP reconciliation. Source hashes are recorded.\\n')
bundle=OUT/'source_bundle.tar.gz'
with tarfile.open(bundle,'w:gz') as tar:
 for p in sorted(OUT.rglob('*')):
  if p.is_file() and p!=bundle:tar.add(p,arcname=str(p.relative_to(OUT)))
print('DELIVERY_DIR='+str(OUT));print('SOURCE_BUNDLE='+str(bundle));print('SOURCE_COUNT='+str(len(sources)))

