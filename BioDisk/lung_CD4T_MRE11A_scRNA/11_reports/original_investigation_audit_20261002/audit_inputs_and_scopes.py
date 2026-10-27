from pathlib import Path
import hashlib,json,re
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,stats
from scipy.io import mmread
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/original_investigation_audit_20261002';R=O/'results';R.mkdir(exist_ok=True)
S1=P/'02_analysis/scenic1_20260929_103513_xRoEnr';M1=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';M2=P/'02_analysis/mechanism_step2_20260927_154919_EsN4gk';T=P/'11_reports/technical_delivery_20261002_ikzf2_stat1_v1/revision_v2/final_release';D=P/'11_reports/stepwise_diagnostic_20261002_v1'
manifest=[]
def rd(p,**kw):
 p=Path(p);manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return pd.read_csv(p,sep='\t',**kw)
def wr(d,n):pd.DataFrame(d).to_csv(R/n,sep='\t',index=False)
mouse_scores=rd(T/'source_tables/M1_sample_scores.tsv');human_scores=rd(T/'source_tables/M2_human_sample_scores.tsv');rna=rd(T/'data/RNA_sample_summary.tsv');end=rd(D/'inputs/receiver_all_positive_cell_samples.tsv')
human_base=rd(M2/'input/human_cells.tsv');hg=rd(M2/'input/human_genes.tsv');hx=mmread(M2/'input/human_counts.mtx').tocsr()
bgset=set('SFTPA1 SFTPA2 SFTPB SFTPC SFTPD SCGB1A1 SCGB3A1 SCGB3A2 ABCA3 SLC34A2 NAPSA LYZ LYZ1 LYZ2 S100A8 S100A9 CSF1R LST1 TYROBP FCER1G'.split())
canon=hg.symbol.astype(str).str.upper();bmask=np.array([v in bgset or bool(re.match(r'^(IGH[VDJMAEG]|IGK[VJC]|IGL[VJC])',v)) for v in canon])
human_base['bg']=np.asarray(hx[bmask].sum(axis=0)).ravel()/human_base.library.to_numpy();human_base['core']=human_base.core.astype(str).str.lower().eq('true')
wr(human_base.groupby(['sample','group','series','subtype','core']).size().reset_index(name='cells'),'human_primary_core_membership.tsv')
del hx
joins=[];human_trace=[];cellchecks=[]
for sp,scopes in [('Mouse',['Core','ActCD4']),('Human',['Core','CD4-IL7R'])]:
 allcells=rd(S1/sp/'All/cells.tsv') if sp=='Mouse' else human_base[human_base.core].copy()
 for arm in ['All','BG1']:
  cells=rd(S1/sp/arm/'cells.tsv')
  for scope in scopes:
   expected=allcells if scope=='Core' else allcells[allcells.subtype==scope]
   if arm=='BG1':expected=expected[expected.bg<=.01]
   z=cells if scope=='Core' else cells[cells.subtype==scope]
   cnt=z.groupby('sample').size();ec=expected.groupby('sample').size()
   qscope=('All' if scope=='Core' else scope) if sp=='Mouse' else ('Core' if scope=='Core' else 'Subtype:'+scope)
   q=mouse_scores[(mouse_scores.dataset=='GSE314943')&(mouse_scores.scope==qscope)&(mouse_scores['filter']==arm)] if sp=='Mouse' else human_scores[(human_scores.scope==qscope)&(human_scores['filter']==arm)]
   for tf in ['IKZF2','STAT1']:
    rr=rna[(rna.species==sp)&(rna.scope==scope)&(rna['filter']==arm)&(rna.tf==tf)]
    final=end[(end.species==sp)&(end.scope==scope)&(end['filter']==arm)&(end.tf==tf)&(end.seed==17)]
    for sample in sorted(set(ec.index)|set(cnt.index)|set(q['sample'])|set(rr['sample'])|set(final['sample'])):
     def val(a,col):
      v=a[a['sample']==sample];return v.iloc[0][col] if len(v) else np.nan
     status='included' if sample in set(final['sample']) else 'upstream_frozen_SCENIC_absent' if sample not in set(cnt.index) else 'module_merge_absent' if sample not in set(q['sample']) else 'RNA_merge_absent'
     joins.append(dict(species=sp,scope=scope,filter=arm,tf=tf,sample=sample,positive_cells_before_frozen_sample_gate=int(ec.get(sample,0)),S1_cells=int(cnt.get(sample,0)),module_cells=val(q,'cells'),RNA_cells=val(rr,'cells'),final_cells=val(final,'receiver_cells'),first_absence=status))
   cellchecks.append(dict(species=sp,scope=scope,filter=arm,pre_sample_gate_cells=len(expected),S1_cells=len(z),frozen_gate_lost_cells=len(expected)-len(z),pre_sample_gate_samples=expected['sample'].nunique(),S1_samples=z['sample'].nunique()))
  if sp=='Human':
   for sample,g in allcells.groupby('sample'):
    bg=g[g.bg<=.01];got=cells[cells['sample']==sample]
    human_trace.append(dict(filter=arm,sample=sample,group=g.group.iloc[0],series=g.series.iloc[0],Core_All=len(g),Core_after_BG1=len(bg),IL7R_All=int((g.subtype=='CD4-IL7R').sum()),IL7R_after_BG1=int((bg.subtype=='CD4-IL7R').sum()),S1_cells=len(got),S1_sample_included=len(got)>0))
wr(joins,'all_merge_anti_join_audit.tsv');wr(cellchecks,'scope_cell_roster_audit.tsv');wr(human_trace,'human_upstream_frozen_sample_gate.tsv')
assert not any(x['first_absence'] in ['module_merge_absent','RNA_merge_absent'] for x in joins)
# Compare original count-layer cell IDs, metadata, raw pseudobulk and later canonical counts.
fp=P/'02_analysis/GSE314943/stage5_20260905_222457_a2dApT/objects/CD4_core.h5ad';a=ad.read_h5ad(fp);x=sparse.csr_matrix(a.layers['counts']);ids=a.obs_names.astype(str);sy=a.var['symbol'].astype(str).str.upper().to_numpy();lib=np.asarray(x.sum(axis=1)).ravel()
m=rd(M1/'comm_input/CD4_cells.tsv');g=rd(M1/'comm_input/CD4_genes.tsv');mx=mmread(M1/'comm_input/CD4_counts.mtx.gz').tocsr().T.tocsr()
assert set(ids)==set(m.cell);ix=pd.Index(ids).get_indexer(m.cell);xx=x[ix]
checks={'cells':len(ids),'cell_ids_exact_set':True,'sample_mapping_equal':bool(np.array_equal(a.obs['sample'].astype(str).to_numpy()[ix],m['sample'])),'group_mapping_equal':bool(np.array_equal(a.obs['group'].astype(str).to_numpy()[ix],m['group'])),'subtype_mapping_equal':bool(np.array_equal(a.obs['Subtype'].astype(str).to_numpy()[ix],m['subtype'])),'library_max_difference':float(np.max(np.abs(lib[ix]-m.library)))}
assert all(checks[k] for k in ['sample_mapping_equal','group_mapping_equal','subtype_mapping_equal']);assert checks['library_max_difference']==0
genes=rd(fp.parents[1]/'tables/genes.tsv');samples=['Saline_'+str(i) for i in [1,2,3]]+['Silica_'+str(i) for i in [1,2,3]]
rawpb=[];expr=[];tfchecks=[]
for scope in ['Core','ActCD4','Naive','Treg']:
 take=np.ones(a.n_obs,bool) if scope=='Core' else a.obs.Subtype.astype(str).eq(scope).to_numpy()
 oldpb=rd(fp.parents[1]/('tables/PB_'+scope+'_counts.tsv.gz'))
 assert list(oldpb.iloc[:,0])==list(a.var_names)
 for sample in samples:
  mask=take&a.obs['sample'].astype(str).eq(sample).to_numpy();calc=np.asarray(x[mask].sum(axis=0)).ravel();delta=np.max(abs(calc-oldpb[sample].to_numpy()))
  rawpb.append(dict(scope=scope,sample=sample,cells=int(mask.sum()),library=int(lib[mask].sum()),max_count_difference=float(delta)));assert delta==0
 for gene in ['IL1R2','PENK','TTC39C','FOXP1','TNFRSF9','LAG3','IKZF2','STAT1']:
  vals=np.asarray(x[:,sy==gene].sum(axis=1)).ravel()
  if gene in ['IKZF2','STAT1']:
   jj=g.canonical.eq(gene).to_numpy();later=np.asarray(mx[:,jj].sum(axis=1)).ravel();delta=float(np.max(abs(vals[ix]-later)));tfchecks.append(dict(gene=gene,max_count_difference=delta));assert delta==0
  for sample in samples:
   mask=take&a.obs['sample'].astype(str).eq(sample).to_numpy()
   raw64=vals[mask].astype('float64')*(1e4/lib[mask].astype('float64'))
   roundtrip=np.expm1(np.log1p(raw64.astype('float32'))).astype('float64')
   expr.append(dict(scope=scope,sample=sample,group=sample.split('_')[0],gene=gene,cells=int(mask.sum()),raw_count=int(vals[mask].sum()),detected_cells=int((vals[mask]>0).sum()),pct=float(100*np.mean(vals[mask]>0)),mean_CP10k=float(np.mean(vals[mask]/lib[mask]*1e4)),mean_CP10k_float64=float(np.mean(raw64)),mean_CP10k_float32_roundtrip=float(np.mean(roundtrip))))
wr(rawpb,'original_count_layer_to_pseudobulk_audit.tsv');wr(expr,'original_count_layer_eight_gene_sample_expression.tsv');wr(tfchecks,'original_to_later_TF_raw_count_audit.tsv')
expr=pd.DataFrame(expr);oldexpr=rd(P/'02_analysis/candidate6_verify_20260914_105332_7aEasw/sample_expression.tsv');oldexpr=oldexpr[(oldexpr.dataset=='GSE314943')&(oldexpr.lineage=='CD4')].copy();oldexpr['scope']=oldexpr.subtype.replace({'All':'Core'})
j=expr.merge(oldexpr,on=['scope','sample','gene'],suffixes=('_recalc','_old'),validate='one_to_one');j['CP10k_delta']=j.mean_CP10k-j.CP10k;j['pct_delta']=j.pct_recalc-j.pct_old
assert len(j)==144 and np.max(abs(j.CP10k_delta))<1e-10 and np.max(abs(j.pct_delta))<1e-10;wr(j,'original_six_gene_expression_reproduction.tsv')
newrna=rna[(rna.species=='Mouse')&(rna['filter']=='All')];j2=expr.merge(newrna,left_on=['scope','sample','gene'],right_on=['scope','sample','tf'],validate='one_to_one',suffixes=('_original','_new'));j2['meanCP10k_delta']=j2.mean_CP10k_original-j2.mean_CP10k_new
j2['float64_meanCP10k_delta']=j2.mean_CP10k_float64-j2.mean_CP10k_new
j2['roundtrip_meanCP10k_delta']=j2.mean_CP10k_float32_roundtrip-j2.mean_CP10k_new
print('TF_MEAN_COMPARISON',len(j2),'raw64 delta',float(np.max(abs(j2.float64_meanCP10k_delta))),'documented roundtrip delta',float(np.max(abs(j2.roundtrip_meanCP10k_delta))),flush=True)
assert len(j2)==len(newrna) and np.max(abs(j2.roundtrip_meanCP10k_delta))<1e-10;wr(j2,'original_TF_count_mean_vs_stepwise_RNA.tsv')
checks.update(original_six_gene_expression_rows=len(j),six_gene_meanCP10k_max_difference=float(np.max(abs(j.CP10k_delta))),TF_meanCP10k_rows=len(j2),TF_float32_vs_new_meanCP10k_max_difference=float(np.max(abs(j2.meanCP10k_delta))),TF_float64_vs_new_meanCP10k_max_difference=float(np.max(abs(j2.float64_meanCP10k_delta))),documented_float32_log1p_expm1_roundtrip_delta=float(np.max(abs(j2.roundtrip_meanCP10k_delta))),downstream_inner_join_losses=0)
(R/'input_scope_audit_summary.json').write_text(json.dumps(checks,indent=2));wr(manifest,'read_input_manifest.tsv')
print(json.dumps(checks));print(pd.DataFrame(human_trace).query("sample in ['CWP149','CWP152']").to_string(index=False));print(pd.DataFrame(joins).query("first_absence != 'included'").to_string(index=False))
