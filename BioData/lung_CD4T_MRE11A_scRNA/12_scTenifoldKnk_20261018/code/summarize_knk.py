from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy.stats import hypergeom,spearmanr
from statsmodels.stats.multitest import multipletests
O=Path('/mnt/BioData/bio/wangmingyuan/lung_CD4T_MRE11A_scRNA/12_scTenifoldKnk_20261018')
V=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA/11_reports/death_priority_20261018_v3/tables')
T=O/'tables';T.mkdir(exist_ok=True)
assert json.loads((O/'compute_complete.json').read_text())['status']=='success'
jobs=json.loads((O/'jobs.json').read_text());mods=pd.read_csv(O/'module_members.tsv',sep='\t');mods=mods[mods.retained.astype(str).str.lower().eq('true')];sets={k:set(v.gene) for k,v in mods.groupby('module')}
allg=[];runs=[];path=[];aud=[]
for j in jobs:
 d=O/'runs'/j['id'];mark=json.loads((d/'complete.json').read_text());assert mark['status']=='success'
 for gene in j['targets']:
  a=mark['targets'][gene];base={k:j[k] for k in ['id','arm','subtype','state','seed','cells','donors','genes']};base['target']=gene
  aud.append(dict(**base,**{k:v for k,v in a.items() if k not in ['gene','genes']},warnings=';'.join(mark['warnings']),seconds=mark['elapsed_seconds']))
  if a['status']!='success':continue
  z=pd.read_csv(d/(gene+'.tsv'),sep='\t');assert z.gene.is_unique and np.isfinite(z.distance).all()
  z=z[z.gene.ne(gene)].copy();univ=set(z.gene);sig=set(z.loc[z['p.value']<.05,'gene']);z['rank_pct']=z.distance.rank(pct=True);z['signed_rank']=z.rank_pct*np.sign(z.directionScore)
  for k,v in base.items():z[k]=v
  allg.append(z);runs.append(dict(**base,nominal=len(sig),fdr=int((z['p.adj']<.05).sum()),up=int((z.direction=='up').sum()),down=int((z.direction=='down').sum()),outgoing=a['outgoing_edges']))
  pp=[]
  for module,members in sets.items():
   ss=(members-{gene})&univ;k=len(ss&sig);n=len(ss);N=len(univ);K=len(sig)
   p=float(hypergeom.sf(k-1,N,n,K)) if n>=3 else np.nan
   pp.append(dict(**base,module=module,defined_members=len(members-{gene}),tested_members=n,hits=k,selected=K,universe=N,fold_enrichment=(k/n)/(K/N) if n and K else np.nan,p=p,hit_genes=';'.join(sorted(ss&sig)),tested_genes=';'.join(sorted(ss)),family=j['id']+'|'+gene+'|module_ORA'))
  valid=[i for i,r in enumerate(pp) if np.isfinite(r['p'])];fdr=multipletests([pp[i]['p'] for i in valid],method='fdr_bh')[1] if valid else []
  for i,q in zip(valid,fdr):pp[i]['FDR']=q
  path.extend(pp)
g=pd.concat(allg,ignore_index=True);r=pd.DataFrame(runs);p=pd.DataFrame(path)
g.to_csv(T/'all_gene_model_results.tsv.gz',sep='\t',index=False);r.to_csv(T/'run_summary.tsv',sep='\t',index=False);p.to_csv(T/'module_ORA.tsv',sep='\t',index=False);pd.DataFrame(aud).to_csv(T/'run_audit.tsv',sep='\t',index=False)
key=['arm','subtype','state','target'];st=[];cons=[]
for kk,z in g.groupby(key):
 base=dict(zip(key,kk));ss=sorted(z.seed.unique());
 for s1,s2 in itertools.combinations(ss,2):
  a=z[z.seed==s1].set_index('gene');b=z[z.seed==s2].set_index('gene');ii=a.index.intersection(b.index);ta=set(a.nlargest(100,'distance').index);tb=set(b.nlargest(100,'distance').index);aa=set(a.index[a['p.value']<.05]);bb=set(b.index[b['p.value']<.05]);st.append(dict(**base,seed1=s1,seed2=s2,distance_spearman=float(spearmanr(a.loc[ii,'distance'],b.loc[ii,'distance']).statistic),top100_jaccard=len(ta&tb)/len(ta|tb),nominal_jaccard=len(aa&bb)/len(aa|bb) if aa|bb else np.nan,direction_agreement=float((a.loc[ii,'direction']==b.loc[ii,'direction']).mean())))
 for gene,zz in z.groupby('gene'):
  cons.append(dict(**base,gene=gene,seeds=len(zz),nominal_seeds=int((zz['p.value']<.05).sum()),fdr_seeds=int((zz['p.adj']<.05).sum()),median_rank=float(zz.rank_pct.median()),median_signed_rank=float(zz.signed_rank.median()),up_seeds=int((zz.direction=='up').sum()),down_seeds=int((zz.direction=='down').sum()),min_model_p=float(zz['p.value'].min()),max_model_p=float(zz['p.value'].max())))
c=pd.DataFrame(cons);c.to_csv(T/'gene_seed_consensus.tsv.gz',sep='\t',index=False);pd.DataFrame(st).to_csv(T/'seed_stability.tsv',sep='\t',index=False)
mp=[]
for kk,z in p.groupby(key+['module']):mp.append(dict(zip(key+['module'],kk))|dict(seeds=len(z),nominal_seeds=int((z.p<.05).sum()),fdr_seeds=int((z.FDR<.05).sum()),members_min=int(z.tested_members.min()),hits_min=int(z.hits.min()),hits_max=int(z.hits.max()),fold_min=float(z.fold_enrichment.min()),fold_max=float(z.fold_enrichment.max()),p_min=float(z.p.min()),p_max=float(z.p.max()),fdr_min=float(z.FDR.min()),fdr_max=float(z.FDR.max())))
pd.DataFrame(mp).to_csv(T/'module_seed_consensus.tsv',sep='\t',index=False)
external=[];eg=[]
for gene,contrasts in {'HIF1A':['GSE253090_Hypoxia','GSE253090_Normoxia','GSE324038_HumanPanT_Hypoxia'],'IL2RA':['GSE121883_spleen','GSE121883_thymus']}.items():
 for contrast in contrasts:
  mouse=contrast.startswith('GSE253090');file=V/(contrast+('_author_DE_mapped.tsv' if mouse else '_DE.tsv'));e=pd.read_csv(file,sep='\t');e=e.dropna(subset=['gene'])
  if mouse:e=e.sort_values(['baseMean','ems'],ascending=[False,True]).drop_duplicates('gene');effect='log2FoldChange'
  else:effect='logFC'
  if effect not in e.columns:raise ValueError((contrast,e.columns.tolist()))
  e=e.set_index('gene');e=e[np.isfinite(e[effect])]
  for kk,z in c[c.target==gene].groupby(key):
   base=dict(zip(key,kk));z=z.set_index('gene');ii=z.index.intersection(e.index);a=z.loc[ii];b=e.loc[ii];stable=a.nominal_seeds.ge(2)&((a.up_seeds==3)|(a.down_seeds==3));selected=ii[stable]
   external.append(dict(**base,contrast=contrast,shared=len(ii),signed_rank_logFC_spearman=float(spearmanr(a.median_signed_rank,b[effect]).statistic),direction_agreement=float((np.sign(a.median_signed_rank)==np.sign(b[effect])).mean()),stable_nominal_shared=len(selected),stable_direction_agreement=float((np.sign(a.loc[selected,'median_signed_rank'])==np.sign(b.loc[selected,effect])).mean()) if len(selected) else np.nan,interpretation='descriptive cross-context concordance; no biological replicate P'))
   for nm in ['BCL2','CASP3','CASP8','CASP9','BAX','BAK1','GSDMD','GSDME','CASP1','PYCARD','NLRP3','LDHA','PKM']:
    if nm in ii:eg.append(dict(**base,contrast=contrast,gene=nm,real_logFC=float(e.loc[nm,effect]),model_nominal_seeds=int(z.loc[nm,'nominal_seeds']),model_up_seeds=int(z.loc[nm,'up_seeds']),model_down_seeds=int(z.loc[nm,'down_seeds'])))
pd.DataFrame(external).to_csv(T/'real_KO_concordance.tsv',sep='\t',index=False);pd.DataFrame(eg).to_csv(T/'key_gene_real_KO_comparison.tsv',sep='\t',index=False)
qa=dict(status='success',networks=len(jobs),knockouts=len(r),gene_rows=len(g),module_tests=int(p.p.notna().sum()),gene_p_family='each single KO/run, package BH',module_fdr_family='eligible prespecified modules per single KO/run',no_biological_P_from_cells=True)
(O/'analysis_complete.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');print(qa,flush=True)
print(pd.DataFrame(mp).query('nominal_seeds>=2').to_string(index=False),flush=True)
exec((O/'code/verify_knk.py').read_text(encoding='utf-8'))
