from pathlib import Path
import sys,json,hashlib
import pandas as pd,numpy as np
from scipy.io import mmread
O=Path(sys.argv[1]);V=O.parent/'integrated_macrophage_CD4_20261015_v4';P=O.parents[1];M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF'
# Reuse validated exact sample-level tests without executing old pipeline.
ns={'__name__':'reused_functions'};exec((V/'code/mouse_stats_v4.py').read_text().split("phen=rd(")[0],ns);assoc=ns['assoc'];S=ns['SAMPLES']
G=['STAT1','IRF1','IRF2','GBP2','GBP5','CASP1','CASP4','CASP11','GSDMD','GSDME','PYCARD','NLRP3','AIM2','IL1B','IL18','CASP3','CASP8','RIPK3','MLKL','ZBP1','NOS2'];L=['TNF','IL15','IL18','IFNG','IFNB1','IFNA1','IFNAB','IL27','EBI3','IL6','OSM','IL12A','IL12B'];R=['IFNGR1','IFNGR2','IFNAR1','IFNAR2','IL27RA','IL6ST','IL6R','OSMR','IL12RB1','IL12RB2','IL2RB','IL2RG','IL15RA','TNFRSF1A','TNFRSF1B']
spec={'genes':G,'ligands':L,'receptors':R,'primary':'reuse original five-gene PyroExec; dissect individual genes','new_branches':'post hoc exploratory; do not replace primary module','unit':'6 matched mice','model':'pooled exact rank correlation and condition-adjusted sensitivity','BH':'all genes within scope,seed,model; original upstream complete-family BH retained','no_mediation':'n=6 insufficient for reliable mediation','new_gene_overlap':'new genes outside excluded panels may occur in original regulons; mark non-independent until rescoring'}
(O/'provenance/analysis_spec.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
manifest=[]
def rd(p):
 p=Path(p);manifest.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});return pd.read_csv(p,sep='\t')
def wr(a,n):a.to_csv(O/'tables'/n,sep='\t',index=False,na_rep='NA')
expr=[];cod=[]
for lin in ['CD4','Mac']:
 md=rd(M/'comm_input'/f'{lin}_cells.tsv');gn=rd(M/'comm_input'/f'{lin}_genes.tsv');fp=M/'comm_input'/f'{lin}_counts.mtx.gz';manifest.append({'path':str(fp),'sha256':hashlib.sha256(fp.read_bytes()).hexdigest()});x=mmread(fp).T.tocsr();genes=pd.Index(gn.canonical);query=G+R if lin=='CD4' else L;query=list(dict.fromkeys(query));ix=genes.get_indexer(query)
 scopes=['Core','Naive','ActCD4','Treg'] if lin=='CD4' else ['AllMac','AMlike','C1qMac','Folr2Mac','Spp1Mac']
 for scope in scopes:
  for sample in S:
   mask=md['sample'].eq(sample).to_numpy() & (np.ones(len(md),bool) if scope in ['Core','AllMac'] else md.subtype.eq(scope).to_numpy());n=int(mask.sum())
   for g,j in zip(query,ix):
    a=x[mask,j].toarray().ravel() if j>=0 and n else None
    expr.append(dict(lineage=lin,scope=scope,sample=sample,group='Silica' if sample.startswith('Silica') else 'Saline',gene=g,cells=n,feature_present=j>=0,CP10k=float(np.mean(a*1e4/md.loc[mask,'library'])) if a is not None else np.nan,positive_cells=int((a>0).sum()) if a is not None else np.nan,pct=float((a>0).mean()*100) if a is not None else np.nan,mouse_symbol=gn.iloc[j].mouse_symbol if j>=0 else 'NOT_IN_MATRIX'))
   if lin=='CD4':
    for name,gg in {'GSDMD_CASP1':['GSDMD','CASP1'],'GSDME_CASP3':['GSDME','CASP3'],'IFNG_receptor':['IFNGR1','IFNGR2'],'typeI_IFN_receptor':['IFNAR1','IFNAR2'],'IL27_receptor':['IL27RA','IL6ST']}.items():
     jj=genes.get_indexer(gg);ok=(jj>=0).all();both=(x[mask][:,jj].toarray()>0).all(axis=1) if ok and n else None;cod.append(dict(scope=scope,sample=sample,pair=name,cells=n,positive_cells=int(both.sum()) if both is not None else np.nan,pct=float(both.mean()*100) if both is not None else np.nan))
expr=pd.DataFrame(expr);wr(expr,'gene_sample_expression.tsv');wr(pd.DataFrame(cod),'same_cell_codetection.tsv')
act=rd(V/'03_receiver/TF_all_regulon_sample_activity.tsv.gz');res=[];de=[];oldpan=rd(M/'tables/module_definition.tsv');excluded=set(oldpan.gene)|{'STAT1','IKZF2'}
for scope in ['Core','Naive','ActCD4','Treg']:
 b=expr[(expr.lineage=='CD4')&(expr.scope==scope)&expr.gene.isin(G)];Y=b.pivot(index='sample',columns='gene',values='CP10k').loc[S].dropna(axis=1)
 d=rd(V/'03_receiver'/f'{scope}_DE_all_genes.tsv.gz');d=d[d.canonical.isin(G)].copy();de.append(d)
 for seed in [17,29,43]:
  X=act[(act.scope==scope)&(act.seed==seed)&(act.regulon=='Stat1(+)')].pivot(index='sample',columns='regulon',values='mean_auc').loc[S]
  a=assoc(X,Y,scope,'none',seed,'STAT1_gene');a['target_excluded_from_frozen_regulons']=a.y.isin(excluded);res.append(a)
res=pd.concat(res);wr(res,'STAT1_gene_associations.tsv');wr(pd.concat(de),'gene_disease_DE_original.tsv')
up=rd(V/'04_chains/all_sender_receiver_associations.tsv.gz');up=up[up.y=='AUC::Stat1(+)'];wr(up,'all_ligand_STAT1_original_families.tsv.gz');wr(up[up.x.isin(L)],'mechanistic_ligand_STAT1.tsv')
# No searching for a desired direction: all nominal alternatives retained, with original BH.
a=up[(up.model=='Pooled')&(up.P<.05)].copy();adj=up[up.model=='ConditionAdjusted'][['scope','sender','x','rho','P','BH']].rename(columns={'rho':'condition_rho','P':'condition_P','BH':'condition_BH'});a=a.merge(adj,on=['scope','sender','x'],validate='1:1');wr(a,'all_nominal_upstream_alternatives.tsv')
net=rd(V/'inputs/TF_target_prior.tsv');wr(net[net.tf.isin(['STAT1','IRF1','IRF2'])&net.target.isin(G)],'signed_prior_edges.tsv')
# Same-state, same-readout links; this table is descriptive post hoc conjunction, not chain P.
down=res[(res.seed==17)&(res.model=='Pooled')&(res.P<.05)&res.y.isin(['GSDMD','CASP1','GSDME','CASP3'])];chain=a.merge(down[['scope','y','rho','P','BH','target_excluded_from_frozen_regulons']],on='scope',suffixes=('_up','_down'));wr(chain,'posthoc_two_segment_candidates.tsv')
wr(pd.DataFrame(manifest).drop_duplicates(),'input_manifest.tsv')
print('DOWNSTREAM_NOMINAL');print(res[(res.seed==17)&(res.model=='Pooled')&(res.P<.05)][['scope','y','rho','P','BH','target_excluded_from_frozen_regulons']].to_string(index=False));print('EXECUTOR_DE');print(pd.concat(de).loc[lambda z:z.canonical.isin(['CASP1','GSDMD','GSDME','CASP3','IRF1']),['scope','canonical','logFC','PValue','FDR']].to_string(index=False));print('CHAIN_ROWS',len(chain))
(O/'provenance/complete.json').write_text(json.dumps({'status':'PASS','association_rows':len(res),'chain_rows':len(chain),'raw_data_cells_as_replicates':False}),encoding='utf-8')