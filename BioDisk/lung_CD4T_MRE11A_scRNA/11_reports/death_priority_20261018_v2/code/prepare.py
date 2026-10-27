import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
from pathlib import Path
import json,hashlib,shutil
import numpy as np,pandas as pd
from scipy import sparse
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA')
O=P/'11_reports/death_priority_20261018_v2';B=P/'11_reports/three_cohort_four_axis_20261017_v16';V=P/'11_reports/death_priority_20261018_v31';A=P/'11_reports/append_analysis_20261017_v18'
for d in ['tables','inputs','logs','objects','code']: (O/d).mkdir(exist_ok=True,parents=True)
def write(df,n):df.to_csv(O/'tables'/n,sep='\t',index=False)
p=pd.read_csv(V/'tables/primary_death_candidates.tsv',sep='\t');q=pd.read_csv(V/'tables/secondary_time_omnibus_candidates.tsv',sep='\t')
p['branch']='primary';q=q.merge(p[['gene','scope']],on=['gene','scope'],how='left',indicator=True);q=q[q._merge=='left_only'].drop(columns='_merge');q['branch']='secondary'
c=pd.concat([p,q],ignore_index=True)[['gene','scope','branch']];write(c,'frozen_candidate_pairs.tsv');cand=set(c.gene)
plan={'version':'independent_route_v2','baseline':'death_priority_20261018_v1; computation v31 unchanged','candidate_pairs':c.to_dict('records'),'primary_selection_unchanged':True,'cohort':'GSE314943 n6 + GSE183682 n4; independent mice n10 strata dataset/day','time':'GSE274372 n12 separate; min10cells sensitivity','statistics':'raw P<.05 nominal; separate BH families by analysis/cohort/variant/model/branch; complete tables; no cell-level inference','candidate_association':'rank residual correlation; exact stratified permutation 2880, condition sensitivity36; leave-one-donor; no causal interpretation; reused same samples are not independent validation','scores':'common measured members across3datasets, exclude all7 new-route candidates, mean log1p sample pseudobulk CP10k, full library denominator; pairwise overlap removed for phenotype-damage; baseline v1 scores preserved','TF':'frozen original motif-supported regulons transferred to common measured genes in joint10; remove all candidate and scoring-panel genes and TF itself; AUCell per cell then donor/state means; >=5 retained targets; no newly inferred joint network','enrichment':'joint disease DE already exists and is reused with original full-family FDR; GO_BP/KEGG ORA rawP<.05 up/down; rank signed sqrt(F) fgsea; separate DB/contrast/scope families','communication':'second-cohort macrophage annotation must precede joint sender claims; same-sample ligand/receptor proxies and NicheNet response prediction; no secretion/causality','sensitivity':['BG1 both cohorts','fixed1/3 Core composition','condition-adjusted associations','leave-one-donor','second cohort day-specific effects descriptive','trajectory K/root baseline retained'],'outputs':'internal biological evidence completion; old baseline preserved; decide additional analyses from actual gaps'}
(O/'protocol.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2),encoding='utf-8')
inputs={};manifest=[]
for ds,sc in [('GSE314943','Core'),('GSE183682','Core'),('GSE274372','TLike')]:
 fs=[B/'inputs'/(ds+'__'+sc+s) for s in ['.cellcounts.npz','.cells.tsv.gz','.genes.txt']]
 x=sparse.load_npz(fs[0]);m=pd.read_csv(fs[1],sep='\t');g=pd.read_csv(fs[2],header=None)[0].astype(str).to_numpy();m['dataset']=ds
 inputs[ds]=(x,m,g)
 for f in fs:manifest.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
common=set.intersection(*(set(g) for x,m,g in inputs.values()))
mods=pd.read_csv(A/'tables/module_member_audit.tsv',sep='\t')[['module','gene']]
auto=pd.read_csv(V/'tables/common_module_members.tsv',sep='\t');auto=auto[auto.module=='Autophagy']
mods=pd.concat([mods,auto]).drop_duplicates();mods['retained']=mods.gene.isin(common)&~mods.gene.isin(cand);write(mods,'module_members.tsv')
sets={mod:sorted(set(z.loc[z.retained,'gene'])) for mod,z in mods.groupby('module')}
ph=['TAct','TReg','Inhib','Naive','Cycle'];damage=['MitoROS','MitoISR','ApoExec','PyroExec'];aux=['OXPHOS','Inflam','Mitophagy','Glycolysis','Autophagy']
coverage=[dict(module=k,defined=mods[mods.module==k].gene.nunique(),retained=len(v)) for k,v in sets.items()];write(pd.DataFrame(coverage),'module_coverage.tsv')
assert all(len(v)>=3 for v in sets.values())
rows=[];expr=[];pairmembers=[]
for a in ph:
 for b in damage:
  sa=set(sets[a])-set(sets[b]);sb=set(sets[b])-set(sets[a])
  for mod,genes in [(a,sa),(b,sb)]:
   pairmembers.extend(dict(pair=a+'__'+b,module=mod,gene=g) for g in sorted(genes))
for ds,(x,m,g) in inputs.items():
 ix={v:i for i,v in enumerate(g)}
 for var in ['Raw','BG1']:
  for scope in (['Core','Naive','ActCD4','Treg'] if ds!='GSE274372' else ['TLike']):
   take=np.ones(len(m),bool) if scope in ['Core','TLike'] else m.subtype.eq(scope).to_numpy()
   if var=='BG1':take&=m.background_fraction.to_numpy()<=.01
   for sam,mm in m.loc[take].groupby('sample'):
    ids=mm.index.to_numpy();cnt=np.asarray(x[ids].sum(0)).ravel();lib=mm.library.sum();cp=cnt/lib*1e4;log=np.log1p(cp)
    r=dict(dataset=ds,sample=sam,scope=scope,variant=var,group=mm.group.iloc[0],day=int(mm.day.iloc[0]),cells=len(mm),library=lib)
    for gene in sorted(cand):r['G_'+gene]=log[ix[gene]];r['CP_'+gene]=cp[ix[gene]]
    for mod,genes in sets.items():r['M_'+mod]=float(log[[ix[v] for v in genes]].mean())
    for a in ph:
     for b in damage:
      for mod,other in [(a,b),(b,a)]:
       genes=sorted(set(sets[mod])-set(sets[other]));r['PAIR_'+a+'__'+b+'__'+mod]=float(log[[ix[v] for v in genes]].mean()) if len(genes)>=3 else np.nan
    rows.append(r)
    for gene in sorted(cand):
     vals=x[ids,ix[gene]].toarray().ravel();expr.append({**{k:r[k] for k in ['dataset','sample','scope','variant','group','day','cells']},'gene':gene,'UMI':int(vals.sum()),'positive_cells':int((vals>0).sum()),'multiUMI_cells':int((vals>=2).sum()),'CP10k':cp[ix[gene]]})
d=pd.DataFrame(rows);fixed=[]
for (ds,var,sam),z in d[d.scope.isin(['Naive','ActCD4','Treg'])].groupby(['dataset','variant','sample']):
 if len(z)!=3:continue
 r=z.iloc[0].to_dict();r.update(scope='Core',variant='FixedComp' if var=='Raw' else 'BG1FixedComp',cells=int(z.cells.sum()),library=z.library.sum())
 for col in d:
  if col.startswith(('M_','CP_','PAIR_')):r[col]=float(z[col].mean())
  elif col.startswith('G_'):r[col]=float(np.log1p(z['CP_'+col[2:]].mean()))
 fixed.append(r)
write(pd.concat([d,pd.DataFrame(fixed)],ignore_index=True),'sample_scores.tsv');write(pd.DataFrame(expr),'candidate_expression.tsv');write(pd.DataFrame(pairmembers),'pairwise_disjoint_module_members.tsv')
de=pd.read_csv(B/'tables/DE_all_models.tsv',sep='\t');write(de[de.gene.isin(cand)],'candidate_existing_DE.tsv')
for sc in ['Core','Naive','ActCD4','Treg']:
 shutil.copy2(B/'tables'/f'joint_stratified__{sc}__Silica_vs_Saline.tsv',O/'inputs'/f'joint_DE_{sc}.tsv')
# Reuse audited implementation, not old candidate lists or old scores.
s=(A/'code/analyze.py').read_text(encoding='utf-8')
head="from pathlib import Path\nimport pandas as pd,numpy as np,itertools,math\nfrom scipy import stats\nO=Path("+repr(str(O))+")\npcache={};results=[];effectrows=[];testids=set()\n"
func=s[s.index('def design('):s.index("ph=plan['phenotypes']")]
fin=s[s.index('def finalize('):s.index("r=finalize(results")]
(O/'code/statslib.py').write_text(head+func+fin,encoding='utf-8')
(O/'inputs/manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('PREPARED',len(c),'gene-state pairs',len(cand),'genes',len(rows)+len(fixed),'score rows',flush=True)
