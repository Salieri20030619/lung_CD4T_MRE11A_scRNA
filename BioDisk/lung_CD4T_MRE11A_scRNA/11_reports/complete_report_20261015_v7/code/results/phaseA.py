import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
from pathlib import Path
import json,hashlib,shutil,itertools,datetime,html,platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O=Path(__file__).resolve().parent.parent
for d in ['tables','figures','provenance','inputs']:(O/d).mkdir(exist_ok=True)
M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF'
D=P/'02_analysis/GSE314943/stage6_20260905_231807_L44XdT/tables'
V=P/'02_analysis/candidate6_verify_20260914_105332_7aEasw'
S=P/'02_analysis/CD4_candidate_slim_20260914_103856_fHVkz4'
mods=['PyroExec','Inflam','ApoExec','OXPHOS','MitoISR','MitoROS','Mitophagy']
genes=['IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1']
manifest=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p,role):
 p=Path(p);manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),role=role))
 return pd.read_csv(p,sep='\t')
def wr(df,n):df.to_csv(O/'tables'/n,sep='\t',index=False,encoding='utf-8',na_rep='NA')
def savefig(fig,n):
 fig.savefig(O/'figures'/(n+'.png'),dpi=170,bbox_inches='tight',facecolor='white');fig.savefig(O/'figures'/(n+'.svg'),bbox_inches='tight',facecolor='white');plt.close(fig)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
# Record specification before calculating results.
spec={'stage':'A only; stop before new target selection','discovery':'GSE314943','frozen_DE':'20260914 inventory, original stage6 coefficients/P/BH unchanged','historical_six':'DE screening plus integrated judgement; not an automatic intersection','DE_nominal_threshold':'P<0.05','historical_effect_threshold':'abs(log2FC)>=0.5','phenotype_modules':mods,'sample_unit':'one equal-weight sample/library summary; never cells','primary_cell_gate':20,'descriptive_cell_gate':1,'primary_filter':'All','sensitivity_filter':'BG1 fixed background panel <=1%','test':'two-sided exact permutation of sample labels using absolute difference in means; at least 2 libraries/group','FDR':'BH separately for Core 7 tests and subtype 21 tests within each filter; differential-response family 21 per filter','animal_mapping':'GEO n=3/group recorded by prior audit; exact animal/pool-to-library mapping unresolved; all new tests are conditional library-level exploratory statistics','excluded':'GSE288887 from discovery/selection; no PRJNA1289989 download','next_stages':'real time and pseudotime remain separate, candidate intersection not run'}
(O/'provenance/specification.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
# Read and preserve required audit context without executing old launchers.
audit_names=['original_investigation_audit_20261002','stepwise_diagnostic_20261002_v1','aligned_evidence_20261002_v1','mouse_sender_20261002_v1','review_R1R2_20261007_v9']
reuse=[]
for name in audit_names:
 root=P/'11_reports'/name
 files=list(root.rglob('*'))
 reuse.append({'directory':str(root),'exists':root.exists(),'file_count':sum(f.is_file() for f in files),'use':{'original_investigation_audit_20261002':'verified original edgeR and six-gene expression reproduction','stepwise_diagnostic_20261002_v1':'sample units, model/threshold distinctions; no old launcher execution','aligned_evidence_20261002_v1':'frozen M1 scores; retain original P/BH families in future association reuse','mouse_sender_20261002_v1':'regulatory/communication input inventory for later phase only','review_R1R2_20261007_v9':'current delivered baseline; CWP evidence is not silica-specific discovery'}[name]})
wr(pd.DataFrame(reuse),'required_report_reuse.tsv')
summary=read(S/'comparison_summary.tsv','20260914 frozen screening counts')
oldhash=read(P/'11_reports/original_investigation_audit_20261002/results/original_source_manifest.tsv','20261002 audit source SHA256')
checks=[];de={};six=[];counts=[]
for scope in ['Core','ActCD4','Naive','Treg']:
 for arm in ['Full','NoBG']:
  f=D/f'{scope}_{arm}_DEG.tsv';a=read(f,'historical DE; preserve original statistical output');de[scope,arm]=a
  a['gene_upper']=a.symbol.str.upper()
  n=int((a.PValue<.05).sum());eligible=(a.PValue<.05)&(a.logFC.abs()>=.5)
  orig=summary[(summary.dataset=='GSE314943')&(summary.scope==scope)&(summary.arm==arm)].iloc[0]
  assert len(a)==orig.rows and n==orig.nominal_p05 and int(eligible.sum())==orig.eligible
  match=oldhash[oldhash.path==str(f)]
  if len(match):assert match.iloc[0].sha256==sha(f)
  checks.append({'file':str(f),'counts_match_20260914':True,'sha256_match_20261002':True if len(match) else 'not_in_manifest'})
  counts.append(dict(scope=scope,arm=arm,tested=len(a),nominal_P05=n,nominal_and_absLFC05=int(eligible.sum()),positive=int((eligible&(a.logFC>0)).sum()),negative=int((eligible&(a.logFC<0)).sum()),original_BH05=int((a.FDR<.05).sum())))
  b=a[a.gene_upper.isin(genes)].copy();b['test_status']='tested'
  missing=[dict(gene_upper=g,symbol=g,test_status='not_tested_in_original_table') for g in genes if g not in set(b.gene_upper)]
  if missing:b=pd.concat([b,pd.DataFrame(missing)],ignore_index=True)
  b['scope']=scope;b['arm']=arm;b['FDR_family']=f'original_{scope}_{arm}_all_{len(a)}_genes';six.append(b)
counts=pd.DataFrame(counts);six=pd.concat(six,ignore_index=True)
assert len(six)==48
wr(counts,'historical_screening_counts.tsv');wr(six,'historical_six_DE_all_scopes.tsv')
(O/'provenance/frozen_DE_verification.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
# All plotted point values remain available.
for scope in ['Core','ActCD4']:
 for arm in ['Full','NoBG']:wr(de[scope,arm],f'volcano_input_{scope}_{arm}.tsv')
fig,axs=plt.subplots(2,2,figsize=(13,10))
for ax,(scope,arm) in zip(axs.flat,itertools.product(['Core','ActCD4'],['Full','NoBG'])):
 a=de[scope,arm];y=-np.log10(a.PValue.clip(lower=np.finfo(float).tiny));nom=(a.PValue<.05)&(a.logFC.abs()>=.5)
 colors=np.where(nom,np.where(a.logFC>0,'#c54636','#286da3'),'#c9cdd1')
 ax.scatter(a.logFC,y,c=colors,s=9,alpha=.65,rasterized=True)
 bh=a.FDR<.05;ax.scatter(a.loc[bh,'logFC'],y[bh],facecolors='none',edgecolors='#222222',s=20,lw=.45)
 for k,g in enumerate(genes):
  r=a[a.gene_upper==g].iloc[0];yy=-np.log10(r.PValue)
  ax.scatter(r.logFC,yy,c='#7c3794',s=27,zorder=4)
  ax.annotate(r.symbol,(r.logFC,yy),xytext=(7,8 if k%2 else -12),textcoords='offset points',fontsize=8,arrowprops=dict(arrowstyle='-',lw=.4))
 ax.axhline(-np.log10(.05),c='grey',ls='--',lw=.7)
 for v in [-.5,.5]:ax.axvline(v,c='grey',ls=':',lw=.7)
 ax.set(title=f'{scope} / {arm} | tested {len(a):,}; P<0.05 {int((a.PValue<.05).sum()):,}',xlabel='log2 fold change: Silica - Saline',ylabel='-log10 original P')
fig.suptitle('Historical GSE314943 DE | original edgeR P and BH retained',fontsize=15)
fig.text(.5,.01,'Color: P<0.05 and |log2FC|>=0.5; black ring: original BH<0.05; purple: historical six (not new selection).',ha='center',fontsize=10)
fig.tight_layout(rect=[0,.04,1,.96]);savefig(fig,'historical_Core_ActCD4_volcano')
fig,ax=plt.subplots(figsize=(13,6));ax.axis('off')
for i,(scope,arm) in enumerate(itertools.product(['Core','ActCD4'],['Full','NoBG'])):
 r=counts[(counts.scope==scope)&(counts.arm==arm)].iloc[0];y=.88-i*.18
 for x,label in zip([.02,.35,.68],[f'{scope} / {arm}\nTested: {r.tested:,}',f'Original P < 0.05\n{r.nominal_P05:,} genes',f'Also |log2FC| >= 0.5\n{r.nominal_and_absLFC05:,} genes']):
  ax.text(x,y,label,transform=ax.transAxes,va='center',bbox=dict(boxstyle='round,pad=.6',fc='#edf3f8',ec='#718ca3'),fontsize=12)
 for x in [.23,.56]:ax.annotate('',xy=(x+.09,y),xytext=(x,y),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=1.3))
ax.text(.5,.10,'Historical six: IL1R2, PENK, TNFRSF9, LAG3, TTC39C, FOXP1\nDE screening + integrated judgement; no documented automatic intersection or unique deterministic rule.',ha='center',va='center',transform=ax.transAxes,fontsize=12,bbox=dict(boxstyle='round,pad=.6',fc='#fff2df',ec='#b69055'))
ax.set_title('Reconstructed historical screening counts (20260914 freeze)',fontsize=16,pad=18);savefig(fig,'historical_screening_flow')
# Heatmap of six historical candidates; avoid implying selection by these panels.
fig,ax=plt.subplots(figsize=(12,6));cols=list(itertools.product(['Core','ActCD4','Naive','Treg'],['Full','NoBG']))
mat=np.array([[six[(six.gene_upper==g)&(six.scope==s)&(six.arm==a)].logFC.iloc[0] for s,a in cols] for g in genes])
im=ax.imshow(mat,cmap='RdBu_r',vmin=-np.nanmax(np.abs(mat)),vmax=np.nanmax(np.abs(mat)),aspect='auto')
for i,g in enumerate(genes):
 for j,(s,a) in enumerate(cols):
  r=six[(six.gene_upper==g)&(six.scope==s)&(six.arm==a)].iloc[0]
  ax.text(j,i,(f'{r.logFC:+.2f}\nP {r.PValue:.2g}\nq {r.FDR:.2g}' if np.isfinite(r.logFC) else 'Not tested'),ha='center',va='center',fontsize=8,color='white' if abs(r.logFC)>np.nanmax(np.abs(mat))*.65 else '#111111')
ax.set(xticks=range(8),xticklabels=[f'{s}\n{a}' for s,a in cols],yticks=range(6),yticklabels=genes,title='Historical six: effect, original P and original BH by scope')
fig.colorbar(im,ax=ax,label='Silica - Saline log2FC');fig.tight_layout();savefig(fig,'historical_six_evidence')
expr=read(V/'sample_expression.tsv','verified six-gene sample expression')
expr=expr[(expr.dataset=='GSE314943')&(expr.lineage=='CD4')&expr.gene.isin(genes)].copy();wr(expr,'historical_six_sample_expression.tsv')
fig,axs=plt.subplots(2,6,figsize=(18,7))
for row,scope in enumerate(['All','ActCD4']):
 for col,g in enumerate(genes):
  ax=axs[row,col];a=expr[(expr.subtype==scope)&(expr.gene==g)]
  for j,gr in enumerate(['Saline','Silica']):
   b=a[a.group==gr].sort_values('sample');ax.scatter(j+np.linspace(-.09,.09,len(b)),b.CP10k,c=['#286da3','#c54636'][j]);ax.hlines(b.CP10k.mean(),j-.15,j+.15,color='black')
  ax.set(title=f'{"Core" if scope=="All" else scope}: {g}',xticks=[0,1],xticklabels=['Saline','Silica'],ylabel='Mean cellular CP10k' if col==0 else '')
fig.suptitle('One point per library; expression summaries do not replace edgeR coefficients',fontsize=15);fig.tight_layout(rect=[0,0,1,.94]);savefig(fig,'historical_six_sample_points')
# Frozen phenotype inputs, verified against the original source.
f=P/'11_reports/aligned_evidence_20261002_v1/inputs/M1_sample_scores.tsv';scores=read(f,'audited sample-level module scores')
original=read(M/'tables/sample_scores.tsv','original M1 score source')
key=['dataset','scope','filter','sample'];sc=['score_'+x for x in mods]
pd.testing.assert_frame_equal(scores.sort_values(key).reset_index(drop=True)[key+sc],original.sort_values(key).reset_index(drop=True)[key+sc],check_exact=False,rtol=1e-12,atol=1e-12)
coverage=read(M/'tables/module_coverage.tsv','module coverage');definitions=read(M/'tables/module_definition.tsv','frozen module membership')
wr(coverage[coverage.dataset=='GSE314943'],'module_coverage_mouse.tsv');wr(definitions[definitions.module.isin(mods)],'module_members.tsv')
mouse=scores[scores.dataset=='GSE314943'].copy();assert not mouse.duplicated(['scope','filter','sample']).any()
wr(mouse[['dataset','scope','filter','sample','group','cells','eligible20']+sc+['detect_'+x for x in mods]],'phenotype_donor_inputs.tsv')
wr(scores[['dataset','scope','filter','sample','group','day','cells','eligible20']],'available_cohort_sample_rosters.tsv')
def effect(a,b,test=True):
 a=np.asarray(a,float);b=np.asarray(b,float);d=b.mean()-a.mean();v=np.r_[a,b];nc=len(a)
 loo=[b.mean()-np.delete(a,i).mean() for i in range(len(a)) if len(a)>1]+[np.delete(b,i).mean()-a.mean() for i in range(len(b)) if len(b)>1]
 r=dict(n_control=len(a),n_case=len(b),control_mean=a.mean(),case_mean=b.mean(),difference=d,case_above_control_mean=int((b>a.mean()).sum()),pair_case_higher=float((b[:,None]>a).mean()),loo_same_sign=float(np.mean(np.sign(loo)==np.sign(d))) if loo else np.nan,loo_min=min(loo) if loo else np.nan,loo_max=max(loo) if loo else np.nan,P=np.nan,permutations=0)
 if test and min(len(a),len(b))>=2:
  null=[]
  for ix in itertools.combinations(range(len(v)),nc):
   mask=np.zeros(len(v),bool);mask[list(ix)]=True;null.append(v[~mask].mean()-v[mask].mean())
  r['P']=float(np.mean(np.abs(null)>=abs(d)-1e-12));r['permutations']=len(null)
 return r
def bh(p):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);ok=np.flatnonzero(np.isfinite(p));idx=ok[np.argsort(p[ok])];q[idx]=np.minimum.accumulate((p[idx]*len(idx)/np.arange(1,len(idx)+1))[::-1])[::-1].clip(0,1);return q
rows=[]
for filt in ['All','BG1']:
 for gate in [1,20]:
  for scope in ['All','Naive','ActCD4','Treg']:
   a=mouse[(mouse['filter']==filt)&(mouse.scope==scope)&(mouse.cells>=gate)]
   for mod in mods:
    x=a[a.group=='Saline']['score_'+mod].dropna();y=a[a.group=='Silica']['score_'+mod].dropna()
    if not len(x) or not len(y):continue
    r=effect(x,y,test=gate==20);r.update(scope='Core' if scope=='All' else scope,filter=filt,min_cells=gate,module=mod,unit='library',evidence='exploratory; independence conditional on historical replicate designation',control_ids=';'.join(a.loc[x.index,'sample']),case_ids=';'.join(a.loc[y.index,'sample']),FDR_family=f'{filt}_gate{gate}_'+('Core_7' if scope=='All' else 'subtypes_21'))
    rows.append(r)
res=pd.DataFrame(rows);res['BH']=res.groupby('FDR_family').P.transform(bh);wr(res,'phenotype_group_comparisons.tsv')
# Test differential responses using within-library subtype contrasts.
inter=[];interinputs=[]
for filt in ['All','BG1']:
 a=mouse[(mouse['filter']==filt)&(mouse.cells>=20)&(mouse.scope!='All')]
 for s1,s2 in itertools.combinations(['Naive','ActCD4','Treg'],2):
  u=a[a.scope==s1].set_index('sample');v=a[a.scope==s2].set_index('sample');ix=u.index.intersection(v.index)
  for mod in mods:
   z=pd.DataFrame({'sample':ix,'group':u.loc[ix,'group'].values,'value':u.loc[ix,'score_'+mod].values-v.loc[ix,'score_'+mod].values}).dropna()
   if z.group.nunique()!=2:continue
   r=effect(z[z.group=='Saline'].value,z[z.group=='Silica'].value);r.update(filter=filt,contrast=f'{s1} minus {s2}',module=mod,FDR_family=f'{filt}_subtype_differential_response_21');inter.append(r)
   z['filter']=filt;z['contrast']=f'{s1} minus {s2}';z['module']=mod;interinputs.append(z)
inter=pd.DataFrame(inter);inter['BH']=inter.groupby('FDR_family').P.transform(bh);wr(inter,'subtype_differential_response.tsv');wr(pd.concat(interinputs),'subtype_differential_response_inputs.tsv')
primary=res[(res['filter']=='All')&(res.min_cells==20)]
fig,axs=plt.subplots(4,7,figsize=(20,11))
for i,scope in enumerate(['All','Naive','ActCD4','Treg']):
 for j,mod in enumerate(mods):
  ax=axs[i,j];a=mouse[(mouse['filter']=='All')&(mouse.scope==scope)]
  for k,gr in enumerate(['Saline','Silica']):
   b=a[a.group==gr].sort_values('sample');pos=k+np.linspace(-.10,.10,len(b));ok=b.cells>=20
   ax.scatter(pos[ok],b.loc[ok,'score_'+mod],c=['#286da3','#c54636'][k],s=28)
   ax.scatter(pos[~ok],b.loc[~ok,'score_'+mod],facecolors='none',edgecolors=['#286da3','#c54636'][k],s=38)
   ax.hlines(b.loc[ok,'score_'+mod].mean(),k-.16,k+.16,color='black')
  name='Core' if scope=='All' else scope;r=primary[(primary.scope==name)&(primary.module==mod)].iloc[0]
  ax.set(title=f'{name}: {mod}\ndelta {r.difference:+.3f}; P {r.P:.2g}; q {r.BH:.2g}',xticks=[0,1],xticklabels=['Saline','Silica'],ylabel='Matched-expression score' if j==0 else '')
fig.suptitle('GSE314943 disease-associated changes within each CD4 scope',fontsize=17)
fig.text(.5,.01,'Each point = one library; hollow point = <20 cells (shown, excluded from test). Exact sample-label permutations; separate Core and subtype BH families.',ha='center',fontsize=11)
fig.tight_layout(rect=[0,.04,1,.95]);savefig(fig,'CD4_program_donor_comparison')
fig,ax=plt.subplots(figsize=(12,5))
m=primary.pivot(index='scope',columns='module',values='difference').reindex(index=['Core','Naive','ActCD4','Treg'],columns=mods)
im=ax.imshow(m,cmap='RdBu_r',vmin=-abs(m.to_numpy()).max(),vmax=abs(m.to_numpy()).max(),aspect='auto')
for i,s in enumerate(m.index):
 for j,mod in enumerate(mods):
  r=primary[(primary.scope==s)&(primary.module==mod)].iloc[0]
  ax.text(j,i,f'{r.difference:+.3f}\nLOO {r.loo_same_sign:.0%}',ha='center',va='center',fontsize=9)
ax.set(xticks=range(7),xticklabels=mods,yticks=range(4),yticklabels=m.index,title='Disease change (Silica - Saline), not absolute score ranking')
fig.colorbar(im,ax=ax,label='Difference in equal-library mean scores');fig.tight_layout();savefig(fig,'CD4_program_change_and_consistency')
roles=pd.DataFrame([
['GSE314943','silica mouse','discovery','3 Saline + 3 Silica libraries','Frozen DE and within-sample subtype program scores; animal/pool mapping not independently resolved'],
['GSE274372','silica mouse','real-time support','D0/D7/D28/D56: 3 libraries each before cell gate','TLike annotation; D0 and D28 only 1/3 library has >=20 cells; revisit sample/state coverage before time inference'],
['GSE183682','silica mouse','direction only','D7/D56: 1 library per condition per time','No replicate-based P values; time and pseudotime are distinct'],
['GSE174725','human BALF silicosis/exposed','limited human support','2 exposed + 3 silicosis in legacy roster','Use revised annotation audit before inference; one case only 5 CD4 cells in legacy scope'],
['GSE296682','human patient','localization only','review existing patient-localization outputs','Not case-control validation'],
['GSE288887','human CWP integrated data','excluded from primary discovery/selection','old report includes CWP plus reused controls','Do not relabel CWP as pure SiO2 silicosis; no PRJNA1289989 download']],columns=['dataset','context','role','sample_design','restriction'])
wr(roles,'cohort_roles_and_reuse.tsv')
# Inventory later-phase inputs without reading raw matrices or running further selection.
inv=[]
for dataset in ['GSE274372','GSE183682','GSE174725']:
 root=P/'02_analysis'/dataset
 for f in root.rglob('*'):
  if f.is_file() and (f.suffix in ['.h5ad','.rds'] or 'pseudotime' in f.name.lower() or (f.suffix=='.tsv' and any(k in f.name for k in ['D7_D0','D28_D0','D56_D0','design','sample','pseudotime']))):
   inv.append(dict(dataset=dataset,path=str(f),bytes=f.stat().st_size,action='inventory only; no new target selection'))
for root in [P/'11_reports/GSE174725_annotation_validation_20261005_v1',P/'11_reports/GSE296682_patient_localization_20261005_v1']:
 for f in root.rglob('*'):
  if f.is_file() and f.suffix in ['.tsv','.json']:inv.append(dict(dataset=root.name,path=str(f),bytes=f.stat().st_size,action='review before later-phase reuse'))
wr(pd.DataFrame(inv),'later_phase_input_inventory.tsv')
# Preserve old methods/QA snippets as separate references.
refs=[M/'methods.txt',D/'methods.txt',P/'11_reports/mouse_sender_20261002_v1/inputs/methods.txt',P/'11_reports/aligned_evidence_20261002_v1/provenance/aligned_methods.json',P/'11_reports/stepwise_diagnostic_20261002_v1/provenance/full_reproduction_audit.json',P/'11_reports/review_R1R2_20261007_v9/package/README_zh.md']
for i,f in enumerate(refs):
 manifest.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f),role='read-only provenance reference'));shutil.copyfile(f,O/'inputs'/f'{i:02d}_{f.name}')
wr(pd.DataFrame(manifest).drop_duplicates('path'),'input_manifest_SHA256.tsv')
# Internal English working report; no customer translation, target shortlist or causal claim.
text=['<!doctype html><html lang="en"><meta charset="utf-8"><title>Plan 1 / phase A</title><style>body{font:15px Arial;max-width:1400px;margin:35px auto;line-height:1.5;color:#172334}table{border-collapse:collapse;font-size:12px}td,th{border:1px solid #ccd4dc;padding:6px}img{max-width:100%}h1,h2{color:#173b58}</style><body>',
'<h1>Silicosis CD4 plan 1: phase A working results</h1><p>New independent revision. Historical outputs preserved. Input inventory, historical screening figures and within-subtype disease program comparisons are complete. New target selection has not started.</p>',
'<h2>The prior report combines human CWP evidence with separately analysed silica mouse evidence</h2><p>The 20261007_R1R2 report describes 19 human donors (10 controls, 9 CWP cases), with cases from GSE227136 and GSE288887. No documented SiO2-specific patient filtering was identified in its method description. Mouse data were analysed separately, not pooled into a single cross-species expression model. The prior human main analysis cannot be represented as a pure silicosis discovery cohort. The new cohort roles below are fixed before new selection.</p>',roles.to_html(index=False),
'<h2>Historical six came from screening and integrated judgement</h2><p>The 20260914 inventory is reproduced from unchanged original edgeR tables. Raw P&lt;0.05 is the nominal screen. The historical |log2FC|&gt;=0.5 filter is shown separately. Volcano plots visualize this evidence; they are not an additional selection algorithm. No unique historical automatic intersection is documented. The six genes are shown for traceability, not guaranteed retention. Full/NoBG are sensitivity arms, not independent replications.</p>',counts.to_html(index=False),
'<img src="figures/historical_screening_flow.png"><img src="figures/historical_Core_ActCD4_volcano.png"><img src="figures/historical_six_evidence.png"><img src="figures/historical_six_sample_points.png">',
'<h2>Program changes are evaluated within each subtype</h2><p>Frozen M1 matched-expression scores are reused and checked against the original source. Scores are per-library means of cellular mean log1p(CP10k) signature expression minus matched control expression. These are custom frozen transcript panels, not protein-level death assays. Seven panels cover pyroptosis/inflammasome, apoptosis, mitochondrial metabolism/stress/antioxidant response, and mitophagy-associated autophagy. Mitophagy is not a comprehensive general-autophagy or flux assay. Module membership and coverage are supplied. Existing panel definitions are reused; no GSE288887 samples enter these computations.</p>',
'<p>Primary comparisons require at least 20 cells per library and two libraries per group. All available-library results at a one-cell gate are descriptive only. Treg Saline_3 has 13 cells; the primary Treg comparison is 2 vs 3, while all-six-library description is retained. BG1 uses the old fixed 1% background panel gate. No threshold is tuned to preserve any gene or subtype. All points and actual sample IDs are exported.</p>',
'<p>New two-sided exact sample-label permutation P values use absolute differences in equal-library means. BH families are Core (7 tests) and subtypes (21 tests), separately by filter. Direct differences in subtype responses are tested on within-library subtype contrasts, with a separate 21-test BH family per filter. Three-vs-three two-sided exact permutations have only 20 assignments and minimum P=0.1; this prevents nominal P&lt;0.05 even for complete separation. Effect sizes and leave-one-library-out direction are therefore central. No significant subtype winner is declared from absolute score or untested effect ranking.</p>',
'<p>Biological independence remains conditional on the historical/GEO replicate designation: exact animal/pool-to-library mapping is not resolved here. All new inferential statistics are explicitly exploratory library-level results; they are not claimed as confirmed six-independent-animal tests. Core also contains composition changes, so it is not a fourth independent subtype.</p>',
'<img src="figures/CD4_program_donor_comparison.png"><img src="figures/CD4_program_change_and_consistency.png">',primary[['scope','module','n_control','n_case','difference','loo_same_sign','P','BH']].to_html(index=False,float_format=lambda x:f'{x:.5g}'),
'<h2>Next decision gate: no new target list yet</h2><p>Use sample-level disease DE, state-appropriate pseudotime and separate real-time evidence to define a prospective intersection. Record eligible universes, exact thresholds and counts at each step. Do not treat cells as biological replicates or infer real time from pseudotime. GSE274372 D0/D28 have limited CD4/TLike coverage under the existing gate, which must be resolved before a valid time-supported intersection. GSE183682 supports directions only. Later evidence cards should connect candidate state, disease effect, temporal direction, regulation/communication, human localization and measurable experimental readouts. IKZF2, STAT1 and the six historical genes have no guaranteed priority.</p>',
'<p>All tables include full numerical values. This is an internal phase-A working report, not a revised customer delivery.</p></body></html>']
(O/'phaseA_working_report.html').write_text('\n'.join(text),encoding='utf-8')
for rec in manifest:assert sha(Path(rec['path']))==rec['sha256']
qa={'status':'PASS','frozen_DE_count_checks':8,'six_DE_rows':len(six),'six_originally_tested_rows':int((six.test_status=='tested').sum()),'M1_original_vs_aligned_scores':'matched','phenotype_rows':len(res),'differential_response_rows':len(inter),'primary_tests':len(primary),'primary_nominal_P05':int((primary.P<.05).sum()),'historical_files_unchanged':True,'new_target_selection_run':False,'new_tests_unit':'library; animal/pool mapping unresolved','python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'provenance/compute_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps(qa));print(primary[['scope','module','difference','loo_same_sign','P','BH']].to_string(index=False))