import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import sys,json,hashlib,itertools,shutil
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy import sparse
O=Path(sys.argv[1]);P=O.parents[1];M=P/'02_analysis/mechanism_step1_20260927_152218_s4V5eF';V=P/'11_reports/customer_sequence_20261015_v2';S=P/'02_analysis/scenic2_20261001_092655_ZKNWD6';D=P/'02_analysis/GSE314943/stage6_20260905_231807_L44XdT/tables'
CAND=['IKZF2','STAT1','IL1R2','PENK','TNFRSF9','LAG3','TTC39C','FOXP1','AHNAK','MDFIC','ST6GALNAC3'];manifest=[]
def rec(f,role):
 f=Path(f);manifest.append(dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),role=role));return f
def rd(f,role,**kw):return pd.read_csv(rec(f,role),sep='\t',**kw)
def wr(a,folder,name):a.to_csv(O/folder/name,sep='\t',index=False,encoding='utf-8',na_rep='NA')
phen=rd(V/'01_CD4_death/all_nonempty_sample_program_scores.tsv','v2 donor-level fixed panels');wr(phen,'01_phenotype','sample_program_scores.tsv')
for f in ['01_CD4_death/disease_program_comparisons.tsv','02_mito_autophagy/disease_program_comparisons.tsv']:
 a=rd(V/f,'unchanged v2 phenotype group statistics');wr(a,'01_phenotype',Path(f).parent.name+'_group_tests.tsv')
lr=pd.read_csv(O/'inputs/LR_prior.tsv',sep='\t');ligands=set(pd.read_csv(O/'inputs/prior_ligand_universe.tsv',sep='\t').gene)&set(lr.ligand);receptors=set(lr.receptor)
complexes={'IL18':(['IL18R1','IL18RAP'],'all'),'IL15':(['IL2RB','IL2RG'],'all'),'TNF':(['TNFRSF1A','TNFRSF1B'],'any'),'TGFB1':(['TGFBR1','TGFBR2'],'all'),'IL33':(['IL1RL1','IL1RAP'],'all'),'IFNG':(['IFNGR1','IFNGR2'],'all')}
frames={};counts={};genes={};expr=[];scope_rows=[];complex_rows=[];universe=[]
for lineage in ['CD4','Mac']:
 md=rd(M/'comm_input'/f'{lineage}_cells.tsv','frozen sample/cell map');gn=rd(M/'comm_input'/f'{lineage}_genes.tsv','frozen canonical gene map');x=mmread(rec(M/'comm_input'/f'{lineage}_counts.mtx.gz','frozen canonical raw counts')).T.tocsr();assert x.shape==(len(md),len(gn))
 frames[lineage]=md;counts[lineage]=x;genes[lineage]=pd.Index(gn.canonical)
 scopes=['Core','Naive','ActCD4','Treg'] if lineage=='CD4' else ['AllMac','AMlike','C1qMac','Folr2Mac','Spp1Mac']
 query=sorted((receptors|set(CAND)|{'IL15RA','IL2RB','IL2RG','IL18R1','IL18RAP'}) if lineage=='CD4' else (ligands|{'IL15RA'}))
 ix=genes[lineage].get_indexer(query);present=ix>=0;xp=x[:,ix[present]].tocsr();yn=xp.multiply((1e4/md.library.to_numpy())[:,None]).tocsr()
 for scope in scopes:
  keep=np.ones(len(md),bool) if scope in ['Core','AllMac'] else md.subtype.eq(scope).to_numpy()
  for sample in sorted(md['sample'].unique()):
   mask=keep&md['sample'].eq(sample).to_numpy();n=int(mask.sum());group=md.loc[md['sample']==sample,'group'].iloc[0]
   scope_rows.append(dict(lineage=lineage,scope=scope,sample=sample,group=group,cells=n,full_library_UMI=md.loc[mask,'library'].sum(),mean_library=md.loc[mask,'library'].mean() if n else np.nan))
   vals={k:np.full(len(query),np.nan) for k in ['CP10k','pct','UMI']}
   if n:
    vals['CP10k'][present]=np.asarray(yn[mask].mean(axis=0)).ravel();vals['pct'][present]=np.asarray((xp[mask]>0).mean(axis=0)).ravel()*100;vals['UMI'][present]=np.asarray(xp[mask].sum(axis=0)).ravel()
   for j,g in enumerate(query):expr.append(dict(lineage=lineage,scope=scope,sample=sample,group=group,cells=n,gene=g,feature_present=bool(present[j]),CP10k=vals['CP10k'][j],pct=vals['pct'][j],UMI=vals['UMI'][j]))
   if lineage=='CD4':
    for lig,(components,rule) in complexes.items():
     ii=genes[lineage].get_indexer(components);ok=(ii>=0).all();det=(x[mask][:,ii].toarray()>0) if n and ok else None
     pct=float((det.all(axis=1) if rule=='all' else det.any(axis=1)).mean()*100) if det is not None else np.nan
     complex_rows.append(dict(ligand=lig,lineage=lineage,scope=scope,sample=sample,group=group,components=';'.join(components),rule=rule,cells=n,joint_pct=pct,measured=bool(ok)))
   elif scope in ['AllMac','AMlike']:
    ii=genes[lineage].get_indexer(['IL15','IL15RA']);ok=(ii>=0).all();det=x[mask][:,ii].toarray()>0 if n and ok else None
    complex_rows.append(dict(ligand='IL15_transpresentation',lineage=lineage,scope=scope,sample=sample,group=group,components='IL15;IL15RA',rule='all',cells=n,joint_pct=float(det.all(axis=1).mean()*100) if det is not None else np.nan,measured=bool(ok)))
expr=pd.DataFrame(expr);wr(expr[expr.lineage=='Mac'],'02_sender','all_ligand_sample_expression.tsv.gz');wr(expr[expr.lineage=='CD4'],'03_receiver','receptor_candidate_sample_expression.tsv.gz');wr(pd.DataFrame(scope_rows),'inputs','sample_scope_roster.tsv');wr(pd.DataFrame(complex_rows),'03_receiver','receptor_complex_codetection.tsv')
# Freeze candidate disease DE and keep original complete-gene BH values.
des=[]
for scope in ['Core','Naive','ActCD4','Treg']:
 f=D/f'{scope}_Full_DEG.tsv';a=rd(f,'original disease DE; all six libraries retained');a['canonical']=a.symbol.str.upper().replace({'MRE11A':'MRE11'});a['scope']=scope;a['origin_family']=f'original_{scope}_Full_all_{len(a)}_genes';wr(a,'03_receiver',f'{scope}_DE_all_genes.tsv.gz');des.append(a)
 for g in CAND:
  b=a[a.canonical==g];universe.append(dict(scope=scope,gene=g,DE_tested=len(b)>0,logFC=b.logFC.iloc[0] if len(b) else np.nan,P=b.PValue.iloc[0] if len(b) else np.nan,BH=b.FDR.iloc[0] if len(b) else np.nan))
wr(pd.DataFrame(universe),'03_receiver','candidate_disease_DE.tsv')
for source in ['AllMac','AMlike']:
 f=P/'02_analysis/mechanism_step3_20260928_103057_lIRL38/tables'/f'DE_Mac_{source}_All.tsv';a=rd(f,'unchanged original six-library sender DE');a['canonical']=a['gene' if 'gene' in a else 'symbol'].str.upper();wr(a,'02_sender',f'{source}_DE_all_genes.tsv.gz')
# All frozen regulons, all nonempty state/sample means. No new SCENIC training.
acts=[];md=frames['CD4'];old=rd(S/'tables/TF_sample_activity.tsv.gz','old sample summaries for roster-preserving parity');coverage=rd(S/'tables/regulon_exclusion_coverage.tsv','frozen regulon sizes')
parity=[]
for seed in [17,29,43]:
 f=S/'objects'/f'Mouse_All_{seed}_DropOverlap_auc.csv.gz';a=pd.read_csv(rec(f,'frozen overlap-removed cell AUCell'),index_col=0);assert a.index.is_unique and set(a.index)==set(md.cell);a=a.loc[md.cell]
 for scope in ['Core','Naive','ActCD4','Treg']:
  keep=np.ones(len(md),bool) if scope=='Core' else md.subtype.eq(scope).to_numpy();b=a.loc[md.loc[keep,'cell']].groupby(md.loc[keep].set_index('cell')['sample']).mean()
  for sample,v in b.iterrows():
   n=int((keep&md['sample'].eq(sample).to_numpy()).sum());group=md.loc[md['sample']==sample,'group'].iloc[0]
   for reg,val in v.items():acts.append(dict(seed=seed,scope=scope,sample=sample,group=group,cells=n,regulon=reg,tf=reg.rsplit('(',1)[0].upper(),mean_auc=val,variant='DropOverlap'))
  oldv=old[(old.species=='Mouse')&(old['filter']=='All')&(old.seed==seed)&(old.variant=='DropOverlap')&(old.scope==scope)]
  for r in oldv.itertuples():
   err=abs(b.loc[r.sample,r.regulon]-r.mean_auc);assert err<1e-12
  parity.append(dict(seed=seed,scope=scope,old_rows=len(oldv),new_rows=b.size,old_nonexcluded_means_equal=True))
acts=pd.DataFrame(acts);wr(acts,'03_receiver','TF_all_regulon_sample_activity.tsv.gz');wr(coverage[(coverage.species=='Mouse')&(coverage['filter']=='All')&(coverage.variant=='DropOverlap')],'03_receiver','TF_regulon_coverage.tsv');wr(pd.DataFrame(parity),'provenance','old_AUC_mean_parity.tsv')
# Eligible LR source pools use detection in >=2 case libraries as a coverage flag, not cell count cutoff.
mac=expr[expr.lineage=='Mac'];recp=expr[expr.lineage=='CD4'];lrrows=[];ligscope=[]
for source in ['AMlike','AllMac']:
 sm=mac[(mac.scope==source)&(mac.group=='Silica')].groupby('gene').agg(sender_case_detected=('UMI',lambda x:int((x>0).sum())),sender_case_mean_pct=('pct','mean'))
 for scope in ['Core','Naive','ActCD4','Treg']:
  rm=recp[(recp.scope==scope)&(recp.group=='Silica')].groupby('gene').agg(receiver_case_detected=('UMI',lambda x:int((x>0).sum())),receiver_case_mean_pct=('pct','mean'))
  b=lr.join(sm,on='ligand').join(rm,on='receptor');b['sender']=source;b['scope']=scope;b['expression_supported']=(b.sender_case_detected>=2)&(b.receiver_case_detected>=2);lrrows.append(b)
  for lig in sorted(ligands):
   v=b[b.ligand==lig];ligscope.append(dict(sender=source,scope=scope,ligand=lig,expression_supported=bool(v.expression_supported.any()),supported_receptors=';'.join(sorted(v.loc[v.expression_supported,'receptor'].unique())),case_source_detected=int(sm.loc[lig,'sender_case_detected']) if lig in sm.index else 0))
wr(pd.concat(lrrows),'02_sender','all_prior_LR_expression_coverage.tsv.gz');wr(pd.DataFrame(ligscope),'02_sender','ligand_scope_eligibility.tsv')
wr(rd(P/'02_analysis/mechanism_step6_20260928_163015_VoRQqA/tables/network_canonical.tsv','frozen DoRothEA mouse network'),'inputs','TF_target_prior.tsv')
wr(pd.DataFrame(manifest).drop_duplicates('path'),'provenance','mouse_input_manifest.tsv')
print('PREPARE_MOUSE_PASS');print('AUC_REGULONS',acts.groupby('seed').regulon.nunique().to_dict());print('LIGAND_ELIGIBILITY',pd.DataFrame(ligscope).groupby(['sender','scope']).expression_supported.sum().to_dict());print('ROSTER',pd.DataFrame(scope_rows).groupby(['lineage','scope']).cells.apply(list).to_dict())