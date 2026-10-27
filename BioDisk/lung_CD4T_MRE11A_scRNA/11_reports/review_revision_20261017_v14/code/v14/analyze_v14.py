import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import sys,json,hashlib,shutil,itertools
import numpy as np,pandas as pd
from scipy.stats import rankdata,false_discovery_control
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');V=P/'11_reports/integrated_macrophage_CD4_20261015_v4';B=P/'11_reports/STAT1_pyro_bridge_20261015_v5';BASE=P/'11_reports/methods_complete_20261016_v13';M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';O=P/'11_reports/review_revision_20261017_v14'
for d in ['supplement_v14','provenance','code','review_inputs']:(O/d).mkdir(parents=True,exist_ok=True)
S=['Saline_1','Saline_2','Saline_3','Silica_1','Silica_2','Silica_3']; scopes=['Core','Naive','ActCD4','Treg']
plan={'date':'2026-10-17','purpose':'response to independent review; targeted post hoc exploratory extension','unit':'six matched biological mice','primary_seed':17,'fixed_axes':{'IL15':['STAT5A','STAT5B'],'TNF':['NFKB1','NFKB2','RELA','RELB'],'IL18':['NFKB1','RELA','JUN','JUNB','FOS','FOSB'],'IL27':['STAT1'],'IFNB1':['STAT1'],'IFNG':['STAT1']},'sources':['AMlike','AllMac'],'scopes':scopes,'BH':'all planned ligand-proximal TF pairs within each source/scope/model, constants NA; existing BH unchanged','comparison_sources':'frozen whole-lung Broad labels; CXCL16 and IFNG, descriptive per mouse; no cell-level testing','receptors':'raw same-cell complex detection; dropout zeros not proof of absence','no_new_screen':'no expansion for significance, no mediation, no model substitution','overlap':'reconstruct DropOverlap and tabulate intersections with fixed modules and controls','original_v13_unchanged':True}
(O/'provenance/analysis_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
manifest=[]
def rd(f):
 f=Path(f);manifest.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()});return pd.read_csv(f,sep='\t')
def wr(a,n):a.to_csv(O/'supplement_v14'/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def bh(p):
 a=np.array(p,float);o=np.full(len(a),np.nan);ok=np.isfinite(a)
 if ok.any():o[ok]=false_discovery_control(a[ok])
 return o
perms=np.array(list(itertools.permutations(range(6))));within=np.array([a+b for a in itertools.permutations(range(3)) for b in itertools.permutations(range(3,6))])
def test(x,y,adjust):
 x=rankdata(x);y=rankdata(y)
 for ix in ([np.arange(3),np.arange(3,6)] if adjust else [np.arange(6)]):x[ix]-=x[ix].mean();y[ix]-=y[ix].mean()
 den=np.linalg.norm(x)*np.linalg.norm(y)
 if den<1e-12:return np.nan,np.nan
 r=x@y/den;null=x[within if adjust else perms]@y/den
 return r,float((np.abs(null)>=abs(r)-1e-12).mean())
# Keep complete original test families, audit BH without reducing to focal rows.
checks=[]
files=[('TF_program',V/'03_receiver/all_TF_program_associations.tsv.gz'),('sender_receiver',V/'04_chains/all_sender_receiver_associations.tsv.gz'),('AUC_disease',V/'03_receiver/all_TF_activity_group_effects.tsv'),('STAT1_gene',B/'tables/STAT1_gene_associations.tsv'),('STAT1_fixed',B/'tables/fixed_composition_STAT1_gene.tsv'),('death',BASE/'supplement/S04_death_program_comparisons.tsv'),('mito',BASE/'supplement/S05_mito_autophagy_comparisons.tsv'),('human',BASE/'supplement/S12_human_disease_comparisons.tsv')]
for label,f in files:
 a=rd(f);wr(a,'S21_complete_'+label+('.tsv.gz' if len(a)>10000 else '.tsv'))
 if 'family' not in a:print('NO_FAMILY',label,list(a));continue
 for fam,z in a.groupby('family',dropna=False):
  if not {'P','BH'}<=set(z):continue
  q=bh(z.P);err=np.nanmax(np.abs(q-z.BH)) if np.isfinite(q).any() else 0
  checks.append(dict(table=label,family=fam,rows=len(z),valid_tests=int(z.P.notna().sum()),max_abs_BH_error=err,status='PASS' if err<1e-10 else 'MISMATCH'))
for sc in scopes:
 a=rd(V/f'03_receiver/{sc}_DE_all_genes.tsv.gz');q=bh(a.PValue);checks.append(dict(table='DE_'+sc,family='all_tested_genes_'+sc,rows=len(a),valid_tests=int(a.PValue.notna().sum()),max_abs_BH_error=float(np.nanmax(abs(q-a.FDR))),status='PASS' if np.nanmax(abs(q-a.FDR))<1e-10 else 'MISMATCH'))
 wr(a,'S21_complete_DE_'+sc+'.tsv.gz')
wr(pd.DataFrame(checks),'S21_BH_recalculation_audit.tsv')
# Preserve rosters, inputs and contamination tables in reviewable supplement.
for f in [V/'inputs/sample_scope_roster.tsv',B/'tables/strict_identity_sensitivity_roster.tsv',B/'tables/CD4_executor_identity_audit.tsv',B/'tables/identity_marker_counts.tsv',B/'tables/source_complex_identity_audit.tsv',B/'tables/gene_sample_expression.tsv',B/'tables/fixed_composition_gene_inputs.tsv',V/'06_support/Core_fixed_composition_TF_inputs.tsv',V/'03_receiver/TF_all_regulon_sample_activity.tsv.gz',V/'02_sender/all_ligand_sample_expression.tsv.gz']:
 a=rd(f);wr(a,'S22_'+f.name)
# Targeted proximal response, preserve full-original AUC group BH.
act=rd(V/'03_receiver/TF_all_regulon_sample_activity.tsv.gz');act=act[act.seed==17];send=rd(V/'02_sender/all_ligand_sample_expression.tsv.gz');extra=rd(B/'tables/gene_sample_expression.tsv');extra=extra[extra.lineage=='Mac']
rows=[];inputrows=[]
for source in plan['sources']:
 for sc in scopes:
  for lig,tfs in plan['fixed_axes'].items():
   sx=send[(send.scope==source)&(send.gene==lig)]
   if len(sx)!=6:sx=extra[(extra.scope==source)&(extra.gene==lig)]
   for tf in tfs:
    sy=act[(act.scope==sc)&(act.tf==tf)]
    ok=len(sx)==6 and len(sy)==6
    x=sx.set_index('sample').reindex(S).CP10k.to_numpy() if len(sx) else np.full(6,np.nan)
    y=sy.set_index('sample').reindex(S).mean_auc.to_numpy() if len(sy) else np.full(6,np.nan)
    ok=ok and np.isfinite(x).all() and np.isfinite(y).all()
    for s,xx,yy in zip(S,x,y):inputrows.append(dict(sender=source,scope=sc,ligand=lig,tf=tf,sample=s,ligand_CP10k=xx,TF_AUC=yy))
    for adj in [False,True]:
     r,p=test(x,y,adj) if ok else (np.nan,np.nan)
     rows.append(dict(sender=source,scope=sc,ligand=lig,tf=tf,model='ConditionAdjusted' if adj else 'Pooled',rho=r,P=p,n_samples=6 if ok else np.nan,status='EVALUABLE' if np.isfinite(p) else 'NOT_EVALUABLE_constant_or_missing'))
a=pd.DataFrame(rows);a['family']=a.sender+'_'+a.scope+'_'+a.model+'_planned_proximal_pairs';a['BH']=a.groupby('family').P.transform(bh);wr(a,'S23_targeted_proximal_associations.tsv');wr(pd.DataFrame(inputrows),'S23_targeted_proximal_sample_inputs.tsv')
g=rd(V/'03_receiver/all_TF_activity_group_effects.tsv');g['tf']=g.regulon.str.split('(').str[0].str.upper();wr(g[(g.seed==17)&g.tf.isin(set(sum(plan['fixed_axes'].values(),[])))],'S23_proximal_AUC_disease_original_family.tsv')
# Receptor complexes and sender trans-presentation from raw matrices.
complexes={'IL15_sender':['IL15','IL15RA'],'IL27_sender':['IL27','EBI3'],'IL15_receiver':['IL2RB','IL2RG'],'IL18_receiver':['IL18R1','IL18RAP'],'IL27_receiver':['IL27RA','IL6ST'],'IFNB1_receiver':['IFNAR1','IFNAR2'],'IFNG_receiver':['IFNGR1','IFNGR2'],'CXCL16_receiver':['CXCR6'],'TNF_receiver':['TNFRSF1A','TNFRSF1B']};co=[];ident=[]
for lin in ['CD4','Mac']:
 md=rd(M/f'comm_input/{lin}_cells.tsv');gn=rd(M/f'comm_input/{lin}_genes.tsv');fp=M/f'comm_input/{lin}_counts.mtx.gz';manifest.append({'path':str(fp),'sha256':hashlib.sha256(fp.read_bytes()).hexdigest()});x=mmread(fp).T.tocsr();genes=pd.Index(gn.canonical)
 for sc in (scopes if lin=='CD4' else ['AllMac','AMlike']):
  for s in S:
   mask=md['sample'].eq(s).to_numpy() & (np.ones(len(md),bool) if sc in ['Core','AllMac'] else md.subtype.eq(sc).to_numpy());n=int(mask.sum())
   for name,gg in complexes.items():
    if (name.endswith('sender'))!=(lin=='Mac'):continue
    ix=genes.get_indexer(gg);ok=(ix>=0).all();b=x[mask][:,ix].toarray()>0 if ok and n else None;z=(b.any(axis=1) if name=='TNF_receiver' else b.all(axis=1)) if b is not None else None
    co.append(dict(lineage=lin,scope=sc,sample=s,group=s.split('_')[0],complex=name,genes=';'.join(gg),rule='any' if name=='TNF_receiver' else 'all',cells=n,positive_cells=int(z.sum()) if z is not None else np.nan,pct=100*z.mean() if z is not None else np.nan,measured=ok))
   if lin=='CD4':
    jj=genes.get_indexer(['GSDMD','SFTPC','SCGB1A1']);z=x[mask][:,jj].toarray();ident.append(dict(scope=sc,sample=s,cells=n,GSDMD_positive=int((z[:,0]>0).sum()),GSDMD_one_UMI=int((z[:,0]==1).sum()),GSDMD_total_UMI=float(z[:,0].sum()),SFTPC_positive=int((z[:,1]>0).sum()),SCGB1A1_positive=int((z[:,2]>0).sum()),GSDMD_with_epithelial=int(((z[:,0]>0)&((z[:,1]>0)|(z[:,2]>0))).sum())))
wr(pd.DataFrame(co),'S24_complex_codetection_raw.tsv');wr(pd.DataFrame(ident),'S22_GSDMD_UMI_source_check.tsv')
# Reconstruct exclusion exactly and enumerate residual controls overlap.
edge=rd(P/'02_analysis/scenic1_20260929_103513_xRoEnr/Mouse/All/regulon_edges.tsv');edge['target_canonical']=edge.target.str.upper().replace({'MRE11A':'MRE11','ATP5A1':'ATP5F1A','ATP5B':'ATP5F1B','ATP5C1':'ATP5F1C','ATP5G1':'ATP5MC1','ATP5O':'ATP5PO','PARK2':'PRKN'});edge['tf_canonical']=edge.tf.str.upper()
mods=rd(BASE/'supplement/S20_all_fixed_module_genes.tsv');ctrl=rd(BASE/'supplement/S20_all_fixed_module_controls.tsv.gz');ctrl=ctrl[ctrl.dataset=='GSE314943'];modall=set(mods.gene);excluded=modall|{'IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','MRE11'}
edge=edge[(~edge.target_canonical.isin(excluded))&(edge.target_canonical!=edge.tf_canonical)];cov=rd(V/'03_receiver/TF_regulon_coverage.tsv');ovs=[]
for (seed,reg),e in edge.groupby(['seed','regulon']):
 tf=e.tf_canonical.iloc[0]
 if tf not in ['STAT1','IKZF2']:continue
 target=set(e.target_canonical);expected=cov[(cov.seed==seed)&(cov.regulon==reg)].retained_genes.iloc[0];assert len(target)==expected,(seed,reg,len(target),expected)
 for mod,z in mods.groupby('module'):
  controls=set(ctrl.loc[ctrl.module==mod,'gene']);ov=target&controls;mg=target&set(z.gene)
  ovs.append(dict(seed=seed,tf=tf,regulon_genes=len(target),module=mod,module_overlap=len(mg),module_overlap_genes=';'.join(sorted(mg)),control_overlap=len(ov),control_overlap_genes=';'.join(sorted(ov)),GSDMD_in_regulon='GSDMD' in target))
wr(pd.DataFrame(ovs),'S25_regulon_module_control_overlap.tsv');wr(edge[edge.tf_canonical.isin(['STAT1','IKZF2'])],'S25_reconstructed_focal_regulon_members.tsv')
# Whole-lung competing sources: retained frozen labels; descriptive only.
import anndata as ad
fp=P/'02_analysis/GSE314943/stage3_20260905_214228_VPtd0r/objects/whole_annotated.h5ad';manifest.append({'path':str(fp),'sha256':hashlib.sha256(fp.read_bytes()).hexdigest()});ad0=ad.read_h5ad(fp);symbols=ad0.var.symbol.astype(str).str.upper();src=[]
for gene in ['CXCL16','IFNG']:
 idx=np.flatnonzero(symbols.eq(gene));assert len(idx)>0
 vals=np.asarray(ad0.layers['counts'][:,idx].sum(axis=1)).ravel();lib=ad0.obs.total_counts.to_numpy();df=ad0.obs[['sample','group','Broad']].copy();df['UMI']=vals;df['CP10k']=vals*1e4/lib
 for (s,grp,label),z in df.groupby(['sample','group','Broad'],observed=True):
  src.append(dict(gene=gene,sample=s,group=grp,source=str(label),cells=len(z),positive_cells=int((z.UMI>0).sum()),positive_pct=float(100*(z.UMI>0).mean()),total_UMI=int(z.UMI.sum()),mean_CP10k=float(z.CP10k.mean()),annotation='frozen Broad; not revalidated here'))
wr(pd.DataFrame(src),'S26_whole_lung_competing_sources.tsv')
wr(pd.DataFrame(manifest).drop_duplicates('path'),'S21_input_manifest.tsv')
status={'BH_families':len(checks),'BH_mismatches':[z for z in checks if z['status']!='PASS'],'new_proximal_rows':len(a),'new_proximal_nominal':int((a.P<.05).sum()),'new_proximal_BH':int((a.BH<.05).sum()),'new_cell_level_tests':0,'python':sys.executable}
(O/'provenance/analysis_complete.json').write_text(json.dumps(status,indent=2),encoding='utf-8');print(json.dumps(status));print('PROXIMAL',a[(a.P<.05)&(a.scope=='ActCD4')].to_string(index=False));print('OVERLAP',pd.DataFrame(ovs).query('seed==17')[['tf','module','control_overlap']].to_string(index=False))