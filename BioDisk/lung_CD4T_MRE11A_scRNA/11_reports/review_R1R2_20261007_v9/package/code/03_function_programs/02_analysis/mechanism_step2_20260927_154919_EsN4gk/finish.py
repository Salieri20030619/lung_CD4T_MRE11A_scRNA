import os,json,tarfile
from pathlib import Path
import pandas as pd
import numpy as np
from shared import effect
r=Path(os.environ['ROOT']);o=Path(os.environ['OUT']);prev=Path(os.environ['PREV'])
def wr(d,n):d.to_csv(o/'tables'/n,sep='\t',index=False,na_rep='NA')
# Mouse module membership effects and leave-one-gene-out direction summaries.
c=pd.read_csv(o/'tables/module_gene_contributions.tsv',sep='\t');l=pd.read_csv(o/'tables/module_leave_one_gene.tsv',sep='\t')
rows=[]
for (scope,filt,mod,g),v in c[(c.dataset=='GSE314943')&c.eligible20].groupby(['scope','filter','module','gene'],observed=True):
    d=effect(v.loc[v.group=='Saline','contribution'],v.loc[v.group=='Silica','contribution'])
    rows.append(dict(scope=scope,filter=filt,module=mod,gene=g,**d))
wr(pd.DataFrame(rows),'mouse_gene_effects.tsv')
rows=[]
for (scope,filt,mod,g),v in l[(l.dataset=='GSE314943')&l.eligible20].groupby(['scope','filter','module','removed_gene'],observed=True):
    a=effect(v.loc[v.group=='Saline','full_score'],v.loc[v.group=='Silica','full_score'])
    b=effect(v.loc[v.group=='Saline','reduced_score'],v.loc[v.group=='Silica','reduced_score'])
    rows.append(dict(scope=scope,filter=filt,module=mod,removed_gene=g,baseline_difference=a['difference'],same_direction=np.sign(a['difference'])==np.sign(b['difference']),**b))
wr(pd.DataFrame(rows),'mouse_leave_one_gene_effects.tsv')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
h=pd.read_csv(o/'tables/human_source_models.tsv',sep='\t');q=h[(h.scope=='Core')&(h['filter']=='All')&(h.model=='Source')&h.metric.str.startswith('score_')].copy()
assert (q.status=='FIT').any(),'No estimable human Core source-adjusted module models'
q['module']=q.metric.str.replace('score_','',regex=False)
fig,ax=plt.subplots(figsize=(6,5));q=q.iloc[::-1]
for i,z in enumerate(q.itertuples()):
    if np.isfinite(z.coefficient):ax.errorbar(z.coefficient,i,xerr=[[z.coefficient-z.CI_low],[z.CI_high-z.coefficient]],fmt='o',color='#176B9B',capsize=2)
ax.axvline(0,color='grey',lw=.8);ax.set(yticks=range(len(q)),yticklabels=q.module.tolist(),xlabel='Adj effect');ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(o/'figures/human_modules.png',dpi=300);plt.close(fig)
(o/'captions.txt').write_text('human_modules.png: y=module abbreviations inherited from step1 module_definition.tsv; x=source-adjusted CWP-Control mean-score coefficient; bars=model-based95% confidence intervals, each subject weighted equally, Core scope, >=20 cells, All filter. Models use group+series; confidence intervals describe these exploratory models. No significance stars. Adj effect=adjusted effect.\nCaspGsd=CASP1 and GSDMD co-detected; Pyro3=CASP1/GSDMD/PYCARD co-detected; AltPyro=CASP3/GSDME co-detected. These are transcript co-detection prerequisites, not protein cleavage measurements.\n')
(o/'methods.txt').write_text('Human lung uses the original final Seurat object and subject/source mapping from existing project scripts. Primary=all final CD4; Core excludes CD4state. No external sample added. The source/condition design and subtype cell counts are exported.\nSame12 frozen modules and matched-expression scoring formula as step1, human controls selected within human cohort then fixed across scopes/filters. Targets excluded. Human and mouse absolute scores are not pooled. Full raw library denominators are retained after symbol aggregation.\nPer-subject mean scores/target CP10k/detection. Source model=response~group+series; SourceDep additionally controls log1p(mean library); WithinG227 compares within the reused mixed-condition source. Minimum20 cells. Non-identifiable or insufficient comparisons explicitly flagged. Original edgeR DEG is not refitted here. Linear models of target expression are exploratory mean-CP10k models, not count-model DEG.\nGroup model BH correction is calculated only when each group has>=10 subjects; all nominal values retained for audit. Partial rank associations require>=10 total subjects,>=4 target-positive subjects and both conditions; condition, source and depth controlled; BH across all target-module associations within scope/filter. Small-sample and source-confounding limitations apply.\nFixed-composition weights are original equal-subject proportions and remain fixed across filters. Missing subtype scores produce NA; low-count subtypes flagged. BG1 uses the same predefined panel fraction<=1% rule.\nGene contributions exactly sum to module mean minus fixed control mean. Leave-one-gene-out keeps original control genes; no new significance testing. Mouse correlation leave-one-sample-out is descriptive with five remaining samples; original six-sample correlations are not selected for significance.\nNicheNet package/cache availability recorded only; ligand-target inference follows verified receiver program definition. No network inference claimed in this stage.\n')
required=['human_source_models.tsv','human_partial_correlations.tsv','module_gene_contributions.tsv','pyro_codetection.tsv','mouse_leave_one_gene_effects.tsv']
assert all((o/'tables'/p).is_file() for p in required)
status=json.loads((o/'python_status.json').read_text());status.update(stage='mechanism_step2_human_and_driver_audit',previous=str(prev),output=str(o),pending=['ligand_target_inference','candidate_regulatory_network'])
(o/'status.json').write_text(json.dumps(status,indent=2))
dest=r/'tmp/lung_mechanism_step2_return.tar.gz';tmp=dest.with_name(dest.name+'.part')
with tarfile.open(tmp,'w:gz') as t:
    for folder in ['tables','figures']:
        for p in sorted((o/folder).iterdir()):t.add(p,arcname=str(p.relative_to(o)))
    for p in sorted(o.iterdir()):
        if p.is_file() and p.suffix in ['.R','.py','.txt','.json','.sh']:t.add(p,arcname=p.name)
assert tmp.stat().st_size<25*1024*1024,'Package exceeds25MiB; all results retained on server'
os.replace(tmp,dest)
(r/'tmp/lung_mechanism_step2_latest.txt').write_text(str(o)+'\n')
print('RETURN_PACKAGE='+str(dest));print('RETURN_BYTES='+str(dest.stat().st_size));print('LUNG_MECHANISM_STEP2_PASS',flush=True)

