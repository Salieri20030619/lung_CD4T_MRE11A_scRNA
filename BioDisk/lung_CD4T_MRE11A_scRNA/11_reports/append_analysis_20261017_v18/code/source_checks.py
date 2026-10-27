from pathlib import Path
import pandas as pd,numpy as np,json,scipy.sparse as sp
P=Path('/mnt/BioDisk/wangmingyuan/lung_CD4T_MRE11A_scRNA');O=P/'11_reports/append_analysis_20261017_v18';B=P/'11_reports/three_cohort_four_axis_20261017_v16';plan=json.loads((O/'protocol_frozen.json').read_text());lr=pd.read_csv(O/'tables/frozen_route_LR_coverage.tsv',sep='\t');x=sp.load_npz(B/'inputs/GSE314943__Core.cellcounts.npz');m=pd.read_csv(B/'inputs/GSE314943__Core.cells.tsv.gz',sep='\t');genes=pd.read_csv(B/'inputs/GSE314943__Core.genes.txt',header=None)[0].tolist();ix={v:i for i,v in enumerate(genes)};tpos=np.asarray(x[:,[ix[v] for v in ['CD3D','CD3E','CD3G','TRAC','CD247'] if v in ix]].sum(1)).ravel()>0;rows=[]
for gene,scope,lig,tf in plan['routes']:
 routes=lr[(lr.dataset=='GSE314943')&(lr.scope==scope)&(lr.sender=='AllMac')&(lr.ligand_parts==lig)].drop_duplicates(['interaction','receptor_parts'])
 if not len(routes):rows.append(dict(gene=gene,scope=scope,ligand=lig,status='NO_FROZEN_EXPRESSED_LR'));continue
 for r in routes.itertuples():
  rec=r.receptor_parts.split(';');mapped=[v for v in rec if v in ix]
  for sample,mm in m.groupby('sample'):
   for var in ['Raw','BG1']:
    z=mm if scope=='Core' else mm[mm.subtype==scope]
    if var=='BG1':z=z[z.background_fraction<=.01]
    ids=z.index.to_numpy();vals=x[ids][:,[ix[v] for v in mapped]].toarray();allpos=(vals>0).all(1) if len(mapped)==len(rec) else np.zeros(len(z),bool);positivecounts=(vals>0).sum(0);frac=(vals>0).mean(0) if len(z) else np.full(len(mapped),np.nan)
    rows.append(dict(gene=gene,scope=scope,ligand=lig,interaction=r.interaction,receptor=r.receptor_parts,sample=sample,group=mm.group.iloc[0],variant=var,cells=len(z),required_subunits=len(rec),measured_subunits=len(mapped),all_subunit_positive_cells=int(allpos.sum()),all_subunit_fraction=float(allpos.mean()) if len(z) else np.nan,all_subunit_Tmarker_positive=int((allpos&tpos[ids]).sum()),subunit_fractions=';'.join(f'{v}:{f:.6g}' for v,f in zip(mapped,frac)),status='OK' if len(z) and len(mapped)==len(rec) else 'NOT_EVALUABLE'))
pd.DataFrame(rows).to_csv(O/'tables/receiver_complex_identity_check.tsv',sep='\t',index=False)
w=pd.read_csv(O/'tables/whole_lung_competing_sources.tsv',sep='\t');source=w[(w.group=='Silica')&(w.cells>=10)].groupby(['gene','Broad']).agg(samples=('sample','nunique'),mean_CP10k=('CP10k','mean'),mean_fraction=('fraction_positive','mean'),min_cells=('cells','min')).reset_index();source['rank_all_eligible']=source.groupby('gene').mean_CP10k.rank(ascending=False,method='min');source['rank_three_samples']=np.nan
for gene,z in source[source.samples==3].groupby('gene'):source.loc[z.index,'rank_three_samples']=z.mean_CP10k.rank(ascending=False,method='min')
source.to_csv(O/'tables/source_competition_summary.tsv',sep='\t',index=False)
# Explicit independence checks of chosen score genes and newly scored TF targets.
mods=pd.read_csv(O/'tables/module_member_audit.tsv',sep='\t');mset=set(mods.loc[~mods.removed_candidate,'gene']);cand=set(g for g,s in plan['focus_pairs']);regs=pd.read_csv(O/'tables/clean_regulon_members.tsv',sep='\t');checks=[]
for tf,z in regs.groupby('tf'):
 keep=set(z.loc[z.retained,'canonical']);assert not keep&mset and not keep&cand;checks.append(dict(tf=tf,retained=len(keep),candidate_overlap=len(keep&cand),programme_overlap=len(keep&mset),status='PASS'))
pd.DataFrame(checks).to_csv(O/'tables/score_independence_audit.tsv',sep='\t',index=False);print('COMPLEX',len(rows),'SOURCES',len(source),'TF_CHECKS',len(checks));(O/'logs/source_complete.txt').write_text('success',encoding='utf-8')