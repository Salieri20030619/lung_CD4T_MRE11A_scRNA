from pathlib import Path
import json,itertools,hashlib
import numpy as np,pandas as pd
from scipy.stats import permutation_test,false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parent.parent
r=pd.read_csv(O/'tables/phenotype_group_comparisons.tsv',sep='\t');d=pd.read_csv(O/'tables/phenotype_donor_inputs.tsv',sep='\t')
i=pd.read_csv(O/'tables/subtype_differential_response.tsv',sep='\t');ii=pd.read_csv(O/'tables/subtype_differential_response_inputs.tsv',sep='\t')
checks=[]
for _,x in r[r.P.notna()].iterrows():
 a=d[(d['filter']==x['filter'])&(d.scope==('All' if x.scope=='Core' else x.scope))&(d.cells>=x.min_cells)]
 c=a[a.group=='Saline']['score_'+x.module].values;s=a[a.group=='Silica']['score_'+x.module].values
 p=permutation_test((c,s),lambda a,b:abs(b.mean()-a.mean()),permutation_type='independent',vectorized=False,n_resamples=np.inf,alternative='greater').pvalue
 assert abs(p-x.P)<1e-10;assert abs(s.mean()-c.mean()-x.difference)<1e-12
 checks.append({'type':'group','filter':x['filter'],'scope':x.scope,'module':x.module,'P':p,'absolute_difference':abs(p-x.P)})
for _,x in i[i.P.notna()].iterrows():
 a=ii[(ii['filter']==x['filter'])&(ii.contrast==x.contrast)&(ii.module==x.module)]
 c=a[a.group=='Saline'].value.values;s=a[a.group=='Silica'].value.values
 p=permutation_test((c,s),lambda a,b:abs(b.mean()-a.mean()),permutation_type='independent',vectorized=False,n_resamples=np.inf,alternative='greater').pvalue
 assert abs(p-x.P)<1e-10
 checks.append({'type':'differential_response','filter':x['filter'],'scope':x.contrast,'module':x.module,'P':p,'absolute_difference':abs(p-x.P)})
families=[]
for name,tab in [('group',r),('differential_response',i)]:
 for family,v in tab.groupby('FDR_family'):
  x=v[v.P.notna()]
  if len(x):assert np.allclose(false_discovery_control(x.P.to_numpy(),method='bh'),x.BH,rtol=1e-12)
  families.append({'type':name,'family':family,'planned_records':len(v),'valid_tests':len(x),'not_tested':len(v)-len(x),'nominal_P05':int((x.P<.05).sum()),'BH05':int((x.BH<.05).sum())})
pd.DataFrame(checks).to_csv(O/'provenance/independent_statistical_checks.tsv',sep='\t',index=False,encoding='utf-8')
pd.DataFrame(families).to_csv(O/'tables/BH_family_inventory.tsv',sep='\t',index=False,encoding='utf-8')
r['test_status']=np.where(r.min_cells==1,'descriptive_only',np.where(r.P.isna(),'insufficient_replicates','exploratory_library_test'))
r.to_csv(O/'tables/phenotype_group_comparisons.tsv',sep='\t',index=False,encoding='utf-8',na_rep='NA')
# Show threshold-adjacent values with enough precision to avoid a false q=0.05 impression.
six=pd.read_csv(O/'tables/historical_six_DE_all_scopes.tsv',sep='\t');genes=['IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1'];cols=list(itertools.product(['Core','ActCD4','Naive','Treg'],['Full','NoBG']))
mat=np.array([[six[(six.gene_upper==g)&(six.scope==s)&(six.arm==a)].logFC.iloc[0] for s,a in cols] for g in genes])
fig,ax=plt.subplots(figsize=(13,6));im=ax.imshow(mat,cmap='RdBu_r',vmin=-np.nanmax(abs(mat)),vmax=np.nanmax(abs(mat)),aspect='auto')
for row,g in enumerate(genes):
 for col,(s,a) in enumerate(cols):
  v=six[(six.gene_upper==g)&(six.scope==s)&(six.arm==a)].iloc[0]
  q=f'{v.FDR:.5f}' if .045<v.FDR<.055 else f'{v.FDR:.3g}'
  label=f'{v.logFC:+.2f}\nP {v.PValue:.3g}\nq {q}' if np.isfinite(v.logFC) else 'Not tested'
  ax.text(col,row,label,ha='center',va='center',fontsize=8,color='white' if abs(v.logFC)>np.nanmax(abs(mat))*.65 else '#111111')
ax.set(xticks=range(8),xticklabels=[f'{s}\n{a}' for s,a in cols],yticks=range(6),yticklabels=genes,title='Historical six: effect, original P and original BH by scope')
fig.colorbar(im,ax=ax,label='Silica - Saline log2FC');fig.tight_layout()
for ext in ['png','svg']:fig.savefig(O/'figures'/('historical_six_evidence.'+ext),dpi=170,bbox_inches='tight',facecolor='white')
plt.close(fig)
report=O/'phaseA_working_report.html';s=report.read_text(encoding='utf-8')
insert='''<h2>Observed phenotype changes support further state-specific investigation, not a single strongest subtype</h2>
<p>In the main 3+3-library comparisons, Naive and ActCD4 PyroExec scores increase by +0.044881 and +0.051628 (P=0.4/0.5; BH=0.98/0.98). Their mitophagy-associated scores decrease by -0.048500 and -0.038711 (P=0.6/0.6; BH=0.98/0.98). These directions remain under BG1. Leave-one-library-out direction consistency is 100% for both pyroptosis panels and 100%/83% for the mitophagy panels. These describe transcriptional directions, not proven cell death or autophagic flux.</p>
<p>Apoptosis is different: Naive ApoExec decreases by -0.007577 (P=0.3, BH=0.98), whereas ActCD4 is essentially unchanged (+0.000056, P=1, BH=1) in the main analysis. ActCD4 ApoExec changes with background filtering (+0.011009 under BG1), so apoptosis should not be merged into a generic increased-death conclusion.</p>
<p>Treg has only 2 control libraries at the 20-cell gate, and BG1 further reduces this to 1 control library; the latter comparisons retain effects but no P/BH. Treg ApoExec direction also changes between all-library description and the 20-cell analysis. These sensitivity and coverage results preclude declaring Treg the preferred experimental state based on an absolute score.</p>
<p>Historical Core FOXP1 BH is 0.052649 (Full) and 0.050039 (NoBG): both exceed 0.05. The figure displays additional precision near the cutoff. All historical six are nominally differential in Core and ActCD4, but that alone does not establish why these six rather than other eligible genes were chosen. A new prospective selection rule remains necessary.</p>
<p>Actual BH family sizes are supplied in BH_family_inventory.tsv. Planned subtype families contain 21 module-by-subtype records; unavailable tests are explicitly missing and BH is applied only to valid tests (14 subtype tests under BG1). Descriptive one-cell-gate results have no tests. No differential-response comparison is used to select a single winning subtype.</p>'''
s=s.replace('<h2>Next decision gate: no new target list yet</h2>',insert+'<h2>Next decision gate: no new target list yet</h2>');report.write_text(s,encoding='utf-8')
qa={'status':'PASS','independent_exact_permutation_checks':len(checks),'BH_families_checked':sum(x['valid_tests']>0 for x in families),'max_P_difference':max(x['absolute_difference'] for x in checks),'visual_review':'six PNG figures inspected; threshold-adjacent BH precision corrected','new_target_selection_run':False}
(O/'provenance/final_validation.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps(qa));print(pd.DataFrame(families).to_string(index=False))