from pathlib import Path
import sys,json,hashlib,itertools
import pandas as pd,numpy as np,anndata as ad
from scipy import sparse
from scipy.stats import permutation_test,spearmanr,rankdata,false_discovery_control
O=Path(sys.argv[1]);V=O.parent/'integrated_macrophage_CD4_20261015_v4';P=O.parents[1];ns={};exec((V/'code/mouse_stats_v4.py').read_text().split('phen=rd(')[0],ns);S=ns['SAMPLES'];assoc=ns['assoc'];bh=ns['bh'];rd=lambda n:pd.read_csv(O/'tables'/n,sep='\t');wr=lambda a,n:a.to_csv(O/'tables'/n,sep='\t',index=False,na_rep='NA')
e=rd('gene_sample_expression.tsv');G=json.loads((O/'provenance/analysis_spec.json').read_text())['genes'];w=pd.read_csv(V/'06_support/Core_fixed_composition_weights.tsv',sep='\t').set_index('state').fixed_weight;Y=None
for s,weight in w.items():
 b=e[(e.scope==s)&e.gene.isin(G)].pivot(index='sample',columns='gene',values='CP10k').reindex(S).dropna(axis=1);Y=b*weight if Y is None else Y+b*weight
X=pd.read_csv(V/'06_support/Core_fixed_composition_TF_inputs.tsv',sep='\t').set_index('sample').loc[S,['Stat1(+)']];out=assoc(X,Y,'Core_fixed','none',17,'STAT1_gene_fixed');wr(out,'fixed_composition_STAT1_gene.tsv');wr(Y.reset_index(),'fixed_composition_gene_inputs.tsv')
# Independent SciPy exact permutation verification, all primary gene associations.
a=rd('STAT1_gene_associations.tsv');act=pd.read_csv(V/'03_receiver/TF_all_regulon_sample_activity.tsv.gz',sep='\t');valid=[]
for row in a[a.seed==17].itertuples():
 x=act[(act.scope==row.scope)&(act.seed==17)&(act.regulon=='Stat1(+)')].set_index('sample').loc[S,'mean_auc'].to_numpy();y=e[(e.scope==row.scope)&(e.gene==row.y)].set_index('sample').loc[S,'CP10k'].to_numpy();rx=rankdata(x);ry=rankdata(y)
 if row.model=='ConditionAdjusted':
  for sl in [slice(0,3),slice(3,6)]:rx[sl]-=rx[sl].mean();ry[sl]-=ry[sl].mean()
  perms=[list(i)+list(j) for i in itertools.permutations(range(3)) for j in itertools.permutations(range(3,6))]
 else:rx-=rx.mean();ry-=ry.mean();perms=list(itertools.permutations(range(6)))
 den=np.linalg.norm(rx)*np.linalg.norm(ry)
 if den<1e-12:assert pd.isna(row.P);continue
 observed=np.dot(rx,ry)/den;pv=np.mean([abs(np.dot(rx[list(ix)],ry)/den)>=abs(observed)-1e-12 for ix in perms]);assert abs(pv-row.P)<1e-12;valid.append(pv)
for _,b in a.groupby('family'):
 ok=b.P.notna();assert np.allclose(false_discovery_control(b.loc[ok,'P']),b.loc[ok,'BH'])
# Human reviewed samples: fixed gene list, mean per-cell normalized expression by donor, no cell gate.
H=P/'11_reports/GSE174725_annotation_validation_20261005_v1/results';hf=H/'broad_annotation_review_SERVER_ONLY.h5ad';mf=H/'frozen_annotation_cells_SERVER_ONLY.tsv.gz';h=ad.read_h5ad(hf);md=pd.read_csv(mf,sep='\t',index_col=0).loc[h.obs_names];x=sparse.csr_matrix(h.layers['counts']);gi=pd.Index(h.var_names.str.upper());rows=[]
for scope in ['Core','MemCD4','NaiveCD4','Treg']:
 keep=md.old_CD4.ne('NotInCore') if scope=='Core' else md.old_CD4.eq(scope)
 for sample in sorted(md['sample'].unique()):
  mask=(keep&md['sample'].eq(sample)).to_numpy();n=int(mask.sum());group=md.loc[md['sample']==sample,'group'].iloc[0]
  for g in G:
   j=gi.get_indexer([g])[0];arr=x[mask,j].toarray().ravel() if n and j>=0 else None
   rows.append(dict(scope=scope,sample=sample,group=group,gene=g,cells=n,CP10k=float(np.mean(arr*1e4/md.loc[mask,'total_counts'])) if arr is not None else np.nan,pct=float((arr>0).mean()*100) if arr is not None else np.nan))
human=pd.DataFrame(rows);wr(human,'human_gene_sample_expression.tsv');stats=[]
for (scope,g),b in human.groupby(['scope','gene']):
 c=b[b.group=='Exposed'].CP10k.dropna().to_numpy();s=b[b.group=='Silicosis'].CP10k.dropna().to_numpy();pv=np.nan
 if min(len(c),len(s))>=2 and np.ptp(np.r_[c,s])>0:pv=permutation_test((c,s),lambda c,s:abs(s.mean()-c.mean()),n_resamples=np.inf,alternative='greater').pvalue
 stats.append(dict(scope=scope,gene=g,n_control=len(c),n_case=len(s),control_mean=np.mean(c) if len(c) else np.nan,case_mean=np.mean(s) if len(s) else np.nan,difference=np.mean(s)-np.mean(c) if len(s) and len(c) else np.nan,P=pv,family='human_'+scope+'_all21'))
stats=pd.DataFrame(stats);stats['BH']=stats.groupby('family').P.transform(bh);wr(stats,'human_gene_effects.tsv')
(O/'provenance/human_sources.json').write_text(json.dumps([{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [hf,mf]],indent=2),encoding='utf-8');(O/'provenance/numeric_validation.json').write_text(json.dumps({'status':'PASS','primary_gene_P_checks':len(valid),'BH_families':a.family.nunique(),'method':'independent scalar dot-product full permutations, scipy BH'}),encoding='utf-8')
print('FIXED');print(out[out.y.isin(['GSDMD','CASP1','GSDME'])][['model','y','rho','P','BH']].to_string(index=False));print('HUMAN');print(stats[stats.gene.isin(['GSDMD','CASP1','GSDME','STAT1'])].to_string(index=False));print('PRIMARY');print(a[(a.seed==17)&a.y.isin(['GSDMD','CASP1'])][['scope','model','y','rho','P','BH']].to_string(index=False))